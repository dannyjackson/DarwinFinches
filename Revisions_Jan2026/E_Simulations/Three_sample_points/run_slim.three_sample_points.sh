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

set -euo pipefail


# ============================================================
# Parse arguments
# ============================================================
#
# Required:
#   -s  selection coefficient
#   -r  per-generation population decline rate
#
# Optional:
#   -f  starting selected-allele frequency at t1 [0.10]
#   -m  msprime recapitation model [hudson]
#   -N  diploid Ne at start [10000]
#   -o  generations from t1 to t2 [10]
#   -O  generations from t1 to t3 [2 * t2 offset]
#   -L  segment length [10000]
#   -n  number of replicates [100]
#
# Example:
#   t1 = 10
#   t2 = 20
#   t3 = 30
#
# sbatch run_simulations.sh \
#     -s 0.05 \
#     -r 0.95 \
#     -f 0.10 \
#     -o 10 \
#     -O 20 \
#     -n 100
# ============================================================

NREPS=100

while getopts "s:r:f:m:N:o:O:L:n:" option; do
    case "${option}" in
        s) S=${OPTARG} ;;
        r) DECLINE_RATE=${OPTARG} ;;
        f) F0=${OPTARG} ;;
        m) MSPRIME_MODEL=${OPTARG} ;;
        N) NE=${OPTARG} ;;
        o) OFFSET2=${OPTARG} ;;
        O) OFFSET3=${OPTARG} ;;
        L) SEG_LEN=${OPTARG} ;;
        n) NREPS=${OPTARG} ;;
        *) echo "Invalid option" >&2; exit 1 ;;
    esac
done


# ----------------------------
# Required arguments
# ----------------------------

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"


# ----------------------------
# Defaults
# ----------------------------

F0="${F0:-0.10}"
MSPRIME_MODEL="${MSPRIME_MODEL:-hudson}"
NE="${NE:-10000}"

# Sampling offsets are measured from t1
OFFSET2="${OFFSET2:-10}"

# If t3 is not supplied, make samples evenly spaced:
# e.g. OFFSET2=10 -> OFFSET3=20
OFFSET3="${OFFSET3:-$((2 * OFFSET2))}"

SEG_LEN="${SEG_LEN:-10000}"


# ----------------------------
# Check sampling order
# ----------------------------

if (( OFFSET2 <= 0 )); then
    echo "ERROR: -o / T2_OFFSET must be > 0" >&2
    exit 1
fi

if (( OFFSET3 <= OFFSET2 )); then
    echo "ERROR: -O / T3_OFFSET must be greater than T2_OFFSET" >&2
    exit 1
fi


# ============================================================
# Paths
# ============================================================

SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"

SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/three_sample_points/simulation_replicates.three_sample_points.slim"

# Change this path if you save the revised Python script
# under a different filename.
PY_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/three_sample_points/three_timepoint_delta_pi_daf.py"

PYTHON="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py/bin/python"


# ============================================================
# Output directory tags
# ============================================================

NE_TAG="Ne_${NE/000/}k"

# Include BOTH sampling offsets so different temporal
# sampling schemes cannot write into the same directory.
GEN_TAG="T2_${OFFSET2}_T3_${OFFSET3}"

AF_TAG="AF_${F0//./_}"


OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/three_sample_points/output/${NE_TAG}/${GEN_TAG}/${AF_TAG}"

mkdir -p "${OUTBASE}/logs"


# One combined TSV per decline-rate / selection-coefficient combo
STATS_TSV="${OUTBASE}/${DECLINE_RATE}/selection_${S}/summary_stats.tsv"

mkdir -p "$(dirname "${STATS_TSV}")"


JOB_ID="${SLURM_JOB_ID:-NA}"

echo "JOB=${JOB_ID} NREPS=${NREPS} S=${S} DECLINE_RATE=${DECLINE_RATE} F0=${F0} NE=${NE} T2_OFFSET=${OFFSET2} T3_OFFSET=${OFFSET3} L=${SEG_LEN}"

echo "OUTBASE=${OUTBASE}"


# ============================================================
# Replicate loop
# ============================================================

for REP in $(seq 1 "${NREPS}"); do

    SEED=$((100000 + REP))

    OUTDIR="${OUTBASE}/${DECLINE_RATE}/selection_${S}/replicate_${REP}"

    mkdir -p "${OUTDIR}"

    cd "${OUTDIR}"


    echo "----------------------------------------"
    echo "REP=${REP} SEED=${SEED} PWD=${PWD}"
    echo "[SLiM] starting..."


    # ========================================================
    # Run SLiM
    # ========================================================

    "${SLIM_BIN}" \
        -d run_id="${REP}" \
        -d seed="${SEED}" \
        -d sel_s="${S}" \
        -d decline_rate="${DECLINE_RATE}" \
        -d f0="${F0}" \
        -d Ne="${NE}" \
        -d L="${SEG_LEN}" \
        -d T2_OFFSET="${OFFSET2}" \
        -d T3_OFFSET="${OFFSET3}" \
        "${SLIM_SCRIPT}" \
        > "run${REP}.log" 2>&1


    echo "[SLiM] done."


    # ========================================================
    # Expected SLiM outputs
    # ========================================================

    TREES="simulation_run${REP}.trees"

    T1_IDS="sample_t1_run${REP}.ids.txt"
    T2_IDS="sample_t2_run${REP}.ids.txt"
    T3_IDS="sample_t3_run${REP}.ids.txt"


    if [[ ! -s "${TREES}" || \
          ! -s "${T1_IDS}" || \
          ! -s "${T2_IDS}" || \
          ! -s "${T3_IDS}" ]]; then

        echo "ERROR: Missing expected SLiM outputs in ${OUTDIR}" >&2

        ls -lh

        exit 2
    fi


    # ========================================================
    # Recapitation + neutral mutations + temporal statistics
    # ========================================================

    echo "[PY] recap + delta-pi / delta-AF stats starting..."


    "${PYTHON}" "${PY_SCRIPT}" \
        --trees "${TREES}" \
        --t1_ids "${T1_IDS}" \
        --t2_ids "${T2_IDS}" \
        --t3_ids "${T3_IDS}" \
        --Ne "${NE}" \
        --mu 2.04e-9 \
        --recomb 1e-8 \
        --L "${SEG_LEN}" \
        --model "${MSPRIME_MODEL}" \
        --seed "${SEED}" \
        --suppress_time_warning \
        --rep "${REP}" \
        --sel_s "${S}" \
        --decline_rate "${DECLINE_RATE}" \
        --offset2 "${OFFSET2}" \
        --offset3 "${OFFSET3}" \
        --tsv_out "${STATS_TSV}" \
        --verbose \
        > "recap_stats_run${REP}.log" 2>&1


    echo "[PY] done."


    # ========================================================
    # Cleanup
    # ========================================================

    # The Python script has already calculated all requested
    # statistics, so the original forward-simulation TS is
    # no longer needed.
    rm -f "${TREES}"

done


echo "ALL DONE (JOB=${JOB_ID}, NREPS=${NREPS})"