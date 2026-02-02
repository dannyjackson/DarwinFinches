#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.err
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --mem=8G
#SBATCH --array=0-9
# 1 s-value * 1000 replicates = 1000 tasks => indices 0..999 --array=0-999

set -euo pipefail

# ----------------------------
# Modules / env
# ----------------------------
module load bcftools
module load vcftools
module load htslib

# Parse arguments
while getopts "s:r:" option; do
    case "${option}" in
        s) S=${OPTARG} ;;
        r) DECLINE_RATE=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"

SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"

# ----------------------------
# Parameters
# ----------------------------
SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/simulation_replicates.slim"
OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/output/"
mkdir -p "$OUTBASE"

# N_REPS=1000
N_REPS=10

TASK_ID="${SLURM_ARRAY_TASK_ID}"
REP=$(( SLURM_ARRAY_TASK_ID % N_REPS + 1 ))

# Per-task output directory
OUTDIR="${OUTBASE}/${DECLINE_RATE}/selection_${S}/replicate_${REP}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "SLURM_JOB_ID=$SLURM_JOB_ID  TASK_ID=$TASK_ID  S=$S  REP=$REP"
echo "PWD=$(pwd)"

# ----------------------------
# Run SLiM
# ----------------------------
SEED=$((100000 + REP))

# Pass selection coefficient into SLiM as SEL_S (must be handled in initialize())
"$SLIM_BIN" -d run_id=$REP -d seed=$SEED -d sel_s=$S -d decline_rate=$DECLINE_RATE -d f0=0.10 \
  "$SLIM_SCRIPT" > "run${REP}.log" 2>&1

slim -d run_id=1 -d seed=100 -d sel_s=0.05 -d decline_rate=0.95 -d f0=0.10 seeded_trees_only.slim


micromamba activate recap_py

python recap_fst_tajd.py \
  --trees output/1/selection_0.00/replicate_1/simulation_run1.trees \
  --t1_ids output/1/selection_0.00/replicate_1/sample_t1_run1.ids.txt \
  --t2_ids output/1/selection_0.00/replicate_1/sample_t2_run1.ids.txt \
  --Ne 150000 --mu 2.04e-9 --recomb 1e-8 --L 50000 \
  --model hudson \
  --seed 1 \
  --out_trees simulation_run1.recap.mut.trees \
  --verbose --suppress_time_warning


# ----------------------------
# Compress + index VCFs
# ----------------------------
bgzip -f "sample_t1_run${REP}.vcf"
tabix -f -p vcf "sample_t1_run${REP}.vcf.gz"

bgzip -f "sample_t2_run${REP}.vcf"
tabix -f -p vcf "sample_t2_run${REP}.vcf.gz"

t1="sample_t1_run${REP}.vcf.gz"
t2="sample_t2_run${REP}.vcf.gz"

# ----------------------------
# Tajima's D
# ----------------------------
# Set L to match your SLiM chromosome length (you had 10000 in the script; adjust as needed)
L=50000

vcftools --gzvcf "$t1" --TajimaD "$L" --out "run${REP}_t1" > "run${REP}_t1.tajima.log" 2>&1
vcftools --gzvcf "$t2" --TajimaD "$L" --out "run${REP}_t2" > "run${REP}_t2.tajima.log" 2>&1

# Append to shared summary safely (avoid race conditions)
SUMDIR="${OUTBASE}/${DECLINE_RATE}/selection_${S}"
mkdir -p "$SUMDIR"
TAJ_SUM="${SUMDIR}/tajimaD.txt"
FST_SUM="${SUMDIR}/fst.txt"
LOCKFILE="${SUMDIR}/.summary.lock"

# Write TajimaD lines (NR==2 contains genome-wide statistic when window==L)
{
  awk -v r="run${REP}" 'NR==2 {print r"\tT1\t"$4}' "run${REP}_t1.Tajima.D"
  awk -v r="run${REP}" 'NR==2 {print r"\tT2\t"$4}' "run${REP}_t2.Tajima.D"
} | flock "$LOCKFILE" -c "cat >> '$TAJ_SUM'"

# ----------------------------
# FST (t1 vs t2)
# ----------------------------
bcftools query -l "$t1" | awk -v r="$REP" '{print $1 "\trun"r"_T1_" $1}' > rename_t1.txt
bcftools query -l "$t2" | awk -v r="$REP" '{print $1 "\trun"r"_T2_" $1}' > rename_t2.txt

bcftools reheader -s rename_t1.txt "$t1" > t1_renamed.vcf.gz
bcftools reheader -s rename_t2.txt "$t2" > t2_renamed.vcf.gz

tabix -f -p vcf t1_renamed.vcf.gz
tabix -f -p vcf t2_renamed.vcf.gz

bcftools query -l t1_renamed.vcf.gz > pop_t1.txt
bcftools query -l t2_renamed.vcf.gz > pop_t2.txt

bcftools merge t1_renamed.vcf.gz t2_renamed.vcf.gz -Oz -o merged.vcf.gz
tabix -f -p vcf merged.vcf.gz

outprefix="run${REP}_t1_vs_t2"
vcftools --gzvcf merged.vcf.gz \
  --weir-fst-pop pop_t1.txt \
  --weir-fst-pop pop_t2.txt \
  --out "$outprefix" \
  > "${outprefix}.log" 2>&1

# Append weighted FST estimate safely
awk -v r="run${REP}" '/weighted Fst estimate/ {print r"\t"$NF}' "${outprefix}.log" \
  | flock "$LOCKFILE" -c "cat >> '$FST_SUM'"

echo "Done: S=$S REP=$REP"
