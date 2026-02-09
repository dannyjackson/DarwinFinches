#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%j.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
# SBATCH --array=1-100   # <-- removed

set -euo pipefail

# ----------------------------
# Parse arguments
# ----------------------------
NREPS=100  # default number of replicates (can override with -n)

while getopts "s:r:f:m:N:o:L:n:" option; do
    case "${option}" in
        s) S=${OPTARG} ;;
        r) DECLINE_RATE=${OPTARG} ;;
        f) F0=${OPTARG} ;;
        m) MSPRIME_MODEL=${OPTARG} ;;
        N) NE=${OPTARG} ;;
        o) OFFSET=${OPTARG} ;;
        L) SEG_LEN=${OPTARG} ;;
        n) NREPS=${OPTARG} ;;   # NEW
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

AF="AF_${F0//./_}"

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"
F0="${F0:-0.10}"
MSPRIME_MODEL="${MSPRIME_MODEL:-hudson}"
NE="${NE:-10000}"
OFFSET="${OFFSET:-10}"
SEG_LEN="${SEG_LEN:-10000}"

# ----------------------------
# Paths
# ----------------------------
SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"
SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/all_parameter_space/simulation_replicates.all_parameter_space.slim"

PY_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/all_parameter_space/recap_fst_tajd.all_parameter_space.py"
PYTHON="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py/bin/python"

NE_TAG="Ne_${NE/000/}k"
GEN_TAG="Gen_${OFFSET}"
AF_TAG="AF_${F0//./_}"

OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/all_parameter_space/output/${NE_TAG}/${GEN_TAG}/${AF_TAG}"
mkdir -p "${OUTBASE}/logs"

# One combined TSV per (decline_rate, selection coeff) combo:
STATS_TSV="${OUTBASE}/${DECLINE_RATE}/selection_${S}/summary_stats.tsv"
LOCKFILE="${STATS_TSV}.lock"
mkdir -p "$(dirname "$STATS_TSV")"

echo "JOB=$SLURM_JOB_ID  NREPS=$NREPS  S=$S  DECLINE_RATE=$DECLINE_RATE  F0=$F0  NE=$NE  OFFSET=$OFFSET  L=$SEG_LEN"
echo "OUTBASE=$OUTBASE"

# ----------------------------
# Replicate loop (replaces array)
# ----------------------------
for REP in $(seq 1 "$NREPS"); do
  SEED=$((100000 + REP))

  OUTDIR="${OUTBASE}/${DECLINE_RATE}/selection_${S}/replicate_${REP}"
  mkdir -p "$OUTDIR"
  cd "$OUTDIR"

  echo "----------------------------------------"
  echo "REP=$REP SEED=$SEED PWD=$PWD"
  echo "[SLiM] starting..."

  "${SLIM_BIN}" \
    -d run_id="${REP}" \
    -d seed="${SEED}" \
    -d sel_s="${S}" \
    -d decline_rate="${DECLINE_RATE}" \
    -d f0="${F0}" \
    -d Ne="${NE}" \
    -d L="${SEG_LEN}" \
    -d T2_OFFSET="${OFFSET}" \
    "${SLIM_SCRIPT}" > "run${REP}.log" 2>&1

  echo "[SLiM] done."

  # Expected outputs from SLiM in this directory:
  TREES="simulation_run${REP}.trees"
  T1_IDS="sample_t1_run${REP}.ids.txt"
  T2_IDS="sample_t2_run${REP}.ids.txt"

  if [[ ! -s "${TREES}" || ! -s "${T1_IDS}" || ! -s "${T2_IDS}" ]]; then
    echo "ERROR: Missing expected SLiM outputs in ${OUTDIR}"
    ls -lh
    exit 2
  fi

  echo "[PY] recap+stats starting..."
  "${PYTHON}" "${PY_SCRIPT}" \
    --trees "${TREES}" \
    --t1_ids "${T1_IDS}" \
    --t2_ids "${T2_IDS}" \
    --Ne "${NE}" --mu 2.04e-9 --recomb 1e-8 --L "${SEG_LEN}" \
    --model "${MSPRIME_MODEL}" \
    --seed "${SEED}" \
    --out_trees "simulation_run${REP}.recap.mut.trees" \
    --suppress_time_warning \
    --rep "${REP}" \
    --sel_s "${S}" \
    --decline_rate "${DECLINE_RATE}" \
    --offset "${OFFSET}" \
    --tsv_out "${STATS_TSV}" \
    --verbose \
    > "recap_stats_run${REP}.log" 2>&1
  echo "[PY] done."
done

echo "ALL DONE (JOB=$SLURM_JOB_ID, NREPS=$NREPS)"
