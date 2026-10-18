#!/bin/bash
#SBATCH --account=mcnew
#SBATCH --job-name=thetas_saf
#SBATCH --partition=standard
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_output/output.thetas_saf.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=100gb

# Step 1 of B2: all-sites, autosome-only site allele frequency likelihoods.
# Self-contained: does NOT source params_preprocessing.sh. That file stops at
# line 20 (LEAD=<SET_VALUE> is a bash syntax error), so MINDEPTHIND was never
# set, and it also runs preprocessing_setup.sh (MSMC clone, pip install, cd ~/programs).
#
# Submit together with step 2 (see B2_EffectivePopulationSize.md):
#   jid=$(sbatch --parsable --array=0-5 B2_1_run_thetas_array.sh)
#   sbatch --dependency=aftercorr:${jid} --array=0-5 B2_2_run_thetas_array.sh

set -euo pipefail

# ---------------- settings (same values as the rest of the pipeline) ----------------
ANGSD=~/programs/angsd/angsd
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches
REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna
CHRLEAD=NC_0          # chromosomes (drops unplaced scaffolds)
SEXCHR=NC_044601      # Z
MINDEPTHIND=4
MININD=4
MINQ=30
MINMAPQ=30

WORKDIR=${OUTDIR}/analyses/effective_pop
SAF_DIR=${WORKDIR}/safs
REGIONS=${WORKDIR}/autosomes.rf.txt

POPS=(crapre crapost forpre forpost parpre parpost)
POP=${POPS[$SLURM_ARRAY_TASK_ID]}
THREADS=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "${WORKDIR}/slurm_output" "$SAF_DIR"
cd "$WORKDIR"

# ---------------- checks ----------------
BAMLIST=${OUTDIR}/referencelists/${POP}bams.txt
[[ -s "$BAMLIST" ]] || { echo "ERROR: missing or empty $BAMLIST" >&2; exit 1; }
[[ -s "${REF}.fai" ]] || { echo "ERROR: missing ${REF}.fai" >&2; exit 1; }
while read -r b; do
  [[ -z "$b" ]] && continue
  [[ -s "$b" ]] || { echo "ERROR: BAM in $BAMLIST not found: $b" >&2; exit 1; }
done < "$BAMLIST"

# Autosome list straight from the reference index (names match the BAM headers).
# Written once; the lock stops array tasks from writing it at the same time.
(
  flock -x 9
  if [[ ! -s "$REGIONS" ]]; then
    awk -v lead="$CHRLEAD" -v z="$SEXCHR" 'index($1, lead) == 1 && index($1, z) != 1 {print $1}' \
      "${REF}.fai" > "${REGIONS}.tmp" && mv "${REGIONS}.tmp" "$REGIONS"
  fi
) 9> "${REGIONS}.lock"

echo "POP=$POP  Threads=$THREADS  BAMs=$(grep -c . "$BAMLIST")  autosomes=$(wc -l < "$REGIONS")"
echo "Filters: -setMinDepthInd $MINDEPTHIND -minInd $MININD -minQ $MINQ -minMapQ $MINMAPQ"

# Remove partial output from an earlier failed run
rm -f "${SAF_DIR}/${POP}.saf.idx" "${SAF_DIR}/${POP}.saf.gz" "${SAF_DIR}/${POP}.saf.pos.gz"

"$ANGSD" \
  -bam "$BAMLIST" \
  -rf "$REGIONS" \
  -out "${SAF_DIR}/${POP}" \
  -dosaf 1 \
  -GL 1 \
  -doMajorMinor 1 \
  -doCounts 1 \
  -setMinDepthInd "$MINDEPTHIND" \
  -minInd "$MININD" \
  -minQ "$MINQ" \
  -minMapQ "$MINMAPQ" \
  -anc "$REF" \
  -nThreads "$THREADS"

[[ -s "${SAF_DIR}/${POP}.saf.idx" ]] || { echo "ERROR: angsd finished but ${POP}.saf.idx is missing" >&2; exit 1; }
echo "Done: ${SAF_DIR}/${POP}.saf.idx"
