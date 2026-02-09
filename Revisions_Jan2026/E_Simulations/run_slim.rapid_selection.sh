#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/output/logs/%x_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --array=1-100

set -euo pipefail

# ----------------------------
# Parse arguments
# ----------------------------
while getopts "s:r:f:m:" option; do
    case "${option}" in
        s) S=${OPTARG} ;;
        r) DECLINE_RATE=${OPTARG} ;;
        f) F0=${OPTARG} ;;
        m) MSPRIME_MODEL=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

AF="AF_${F0//./_}"

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"
F0="${F0:-0.10}"
MSPRIME_MODEL="${MSPRIME_MODEL:-hudson}"

# ----------------------------
# Paths
# ----------------------------
SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"
SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/simulation_replicates.rapid_selection.slim"

PY_ENV="recap_py"
PY_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/recap_fst_tajd.py"

OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/rapid_selection/${AF}/output"
mkdir -p "${OUTBASE}/logs"

# One combined TSV per (decline_rate, selection coeff) combo:
STATS_TSV="${OUTBASE}/${DECLINE_RATE}/selection_${S}/summary_stats.tsv"
LOCKFILE="${STATS_TSV}.lock"

# ----------------------------
# Replicate mapping
# ----------------------------
REP="${SLURM_ARRAY_TASK_ID}"   # because array=1..N
SEED=$((100000 + REP))

OUTDIR="${OUTBASE}/${DECLINE_RATE}/selection_${S}/replicate_${REP}"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "JOB=$SLURM_JOB_ID  ARRAY_TASK=$SLURM_ARRAY_TASK_ID  REP=$REP  S=$S  DECLINE_RATE=$DECLINE_RATE  F0=$F0"
echo "PWD=$PWD"

# ----------------------------
# Run SLiM
# ----------------------------
echo "[SLiM] starting..."
"${SLIM_BIN}" \
  -d run_id="${REP}" \
  -d seed="${SEED}" \
  -d sel_s="${S}" \
  -d decline_rate="${DECLINE_RATE}" \
  -d f0="${F0}" \
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

# Ensure stats directory exists
mkdir -p "$(dirname "$STATS_TSV")"

PYTHON="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py/bin/python"

echo "[PY] recap+stats starting..."
"${PYTHON}" "${PY_SCRIPT}" \
  --trees "${TREES}" \
  --t1_ids "${T1_IDS}" \
  --t2_ids "${T2_IDS}" \
  --Ne 10000 --mu 2.04e-9 --recomb 1e-8 --L 10000 \
  --model "${MSPRIME_MODEL}" \
  --seed "${SEED}" \
  --out_trees "simulation_run${REP}.recap.mut.trees" \
  --suppress_time_warning \
  --rep "${REP}" \
  --sel_s "${S}" \
  --decline_rate "${DECLINE_RATE}" \
  --tsv_out "${STATS_TSV}" \
  --verbose \
  > "recap_stats_run${REP}.log" 2>&1
echo "[PY] done."



# ----------------------------
# Make the TSV append safe under parallel array jobs
# (Only needed if your python writes directly without locking.)
# If you implement TSV writing inside python exactly as shown above, it will still
# usually be okay, but locking is safest.
# ----------------------------
# If you want locking at the bash level instead, remove --tsv_out from python and do:
#   last_line=$(tail -n 1 stats_from_python.tsv)
# Here we rely on python writing the TSV. If you want strict locking:
if command -v flock >/dev/null 2>&1; then
  # create lockfile directory if needed
  touch "${LOCKFILE}"
  # nothing to do; python already appended. This is just a placeholder to show locking.
  :
fi

echo "ALL DONE for REP=${REP}"
