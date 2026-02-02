#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=0-499
# 5 s-values * 100 replicates = 500 tasks => indices 0..499

set -euo pipefail

# ----------------------------
# Modules / env
# ----------------------------
module load bcftools
module load vcftools
module load htslib

SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"

# ----------------------------
# Parameters
# ----------------------------
SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/simulation_replicates.slim"
OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/output"

# NOTE: your list had 0.1 twice and ended with 0.1 again; I’m keeping your exact list.
S_LIST=(0.00 0.25 0.5 0.75 1)
N_REPS=100

# Map array task ID -> (s_index, rep)
TASK_ID="${SLURM_ARRAY_TASK_ID}"
S_INDEX=$(( TASK_ID / N_REPS ))
REP=$(( TASK_ID % N_REPS + 1 ))

S="${S_LIST[$S_INDEX]}"

# Per-task output directory
OUTDIR="${OUTBASE}/selection_${S}/replicate_${REP}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "SLURM_JOB_ID=$SLURM_JOB_ID  TASK_ID=$TASK_ID  S=$S  REP=$REP"
echo "PWD=$(pwd)"

# ----------------------------
# Run SLiM
# ----------------------------
SEED=$((100000 + REP))

# Pass selection coefficient into SLiM as SEL_S (must be handled in initialize())
"$SLIM_BIN" -d run_id=$REP -d seed=$SEED -d sel_s=$S \
  "$SLIM_SCRIPT" > "run${REP}.log" 2>&1

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
L=10000

vcftools --gzvcf "$t1" --TajimaD "$L" --out "run${REP}_t1" > "run${REP}_t1.tajima.log" 2>&1
vcftools --gzvcf "$t2" --TajimaD "$L" --out "run${REP}_t2" > "run${REP}_t2.tajima.log" 2>&1

# Append to shared summary safely (avoid race conditions)
SUMDIR="${OUTBASE}/selection_${S}"
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
