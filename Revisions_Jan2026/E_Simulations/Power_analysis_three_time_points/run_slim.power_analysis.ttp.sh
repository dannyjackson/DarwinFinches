#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp/output/logs/%x_%j.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp/output/logs/%x_%j.err
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=2
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G

set -euo pipefail


# ============================================================
# Parse arguments
#
# REQUIRED:
#
#   -s  selection coefficient
#   -r  population decline rate
#   -f  starting selected-allele frequency
#   -m  msprime recapitation model
#   -N  effective population size
#   -o  generations from T1 -> T2
#   -O  generations from T1 -> T3
#   -L  segment length
#   -n  number of replicates
#
# Example:
#
# sbatch run_slim_power_analysis_ttp.ttp.sh \
#   -s 0.05 \
#   -r 0.95 \
#   -f 0.10 \
#   -m hudson \
#   -N 150000 \
#   -o 10 \
#   -O 20 \
#   -L 10000 \
#   -n 100
# ============================================================

S=""
DECLINE_RATE=""
F0=""
MSPRIME_MODEL=""
NE=""
OFFSET2=""
OFFSET3=""
SEG_LEN=""
NREPS=""

while getopts "s:r:f:m:N:o:O:L:n:" option; do
    case "${option}" in
        s) S="${OPTARG}" ;;
        r) DECLINE_RATE="${OPTARG}" ;;
        f) F0="${OPTARG}" ;;
        m) MSPRIME_MODEL="${OPTARG}" ;;
        N) NE="${OPTARG}" ;;
        o) OFFSET2="${OPTARG}" ;;
        O) OFFSET3="${OPTARG}" ;;
        L) SEG_LEN="${OPTARG}" ;;
        n) NREPS="${OPTARG}" ;;
        *)
            echo "ERROR: Invalid option." >&2
            exit 1
            ;;
    esac
done


# ============================================================
# Require ALL parameters
# ============================================================

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"
: "${F0:?ERROR: -f <f0> is required}"
: "${MSPRIME_MODEL:?ERROR: -m <msprime_model> is required}"
: "${NE:?ERROR: -N <Ne> is required}"
: "${OFFSET2:?ERROR: -o <T1_to_T2_offset> is required}"
: "${OFFSET3:?ERROR: -O <T1_to_T3_offset> is required}"
: "${SEG_LEN:?ERROR: -L <segment_length> is required}"
: "${NREPS:?ERROR: -n <number_of_replicates> is required}"


# ============================================================
# Basic validation
# ============================================================

if (( OFFSET2 <= 0 )); then
    echo "ERROR: T1->T2 offset (-o) must be > 0." >&2
    exit 1
fi

if (( OFFSET3 <= OFFSET2 )); then
    echo "ERROR: T1->T3 offset (-O) must be greater than T1->T2 offset (-o)." >&2
    exit 1
fi

if (( NE <= 0 )); then
    echo "ERROR: Ne (-N) must be > 0." >&2
    exit 1
fi

if (( SEG_LEN <= 0 )); then
    echo "ERROR: segment length (-L) must be > 0." >&2
    exit 1
fi

if (( NREPS <= 0 )); then
    echo "ERROR: number of replicates (-n) must be > 0." >&2
    exit 1
fi


# ============================================================
# Paths
# ============================================================

SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"

SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp/simulation_replicates.power_analysis_ttp.slim"

PY_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp/three_timepoint_delta_pi_daf.py"

PYTHON="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py/bin/python"


# ============================================================
# Output tags
# ============================================================

NE_TAG="Ne_${NE/000/}k"

# Both temporal offsets are encoded in directory name.
#
# Example:
#   OFFSET2=10
#   OFFSET3=20
#
# gives:
#   Gen_10_20
#
GEN_TAG="Gen_${OFFSET2}_${OFFSET3}"

AF_TAG="AF_${F0//./_}"


OUTBASE="/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp/${NE_TAG}/${GEN_TAG}"

mkdir -p "${OUTBASE}/logs"


# ============================================================
# Combined statistics TSV
#
# One TSV per:
#   Ne
#   timing combination
#   decline rate
#   selection coefficient
#   starting allele frequency
#   msprime model
# ============================================================

STATS_TSV="${OUTBASE}/d${DECLINE_RATE}_s${S}_af${F0}_model${MSPRIME_MODEL}.summary_stats.tsv"

mkdir -p "$(dirname "${STATS_TSV}")"


# ============================================================
# Job summary
# ============================================================

echo "============================================================"
echo "JOB=${SLURM_JOB_ID}"
echo "NREPS=${NREPS}"
echo "S=${S}"
echo "DECLINE_RATE=${DECLINE_RATE}"
echo "F0=${F0}"
echo "NE=${NE}"
echo "OFFSET2=${OFFSET2}  (T1 -> T2)"
echo "OFFSET3=${OFFSET3}  (T1 -> T3)"
echo "L=${SEG_LEN}"
echo "MSPRIME_MODEL=${MSPRIME_MODEL}"
echo "OUTBASE=${OUTBASE}"
echo "STATS_TSV=${STATS_TSV}"
echo "============================================================"


# ============================================================
# Replicate loop
# ============================================================

for REP in $(seq 1 "${NREPS}"); do

    SEED=$((100000 + REP))


    # --------------------------------------------------------
    # Keep temporary replicate directories isolated by
    # decline, selection, allele frequency, and model.
    #
    # This prevents simultaneously running parameter
    # combinations from writing into the same replicate dirs.
    # --------------------------------------------------------

    RUNBASE="${OUTBASE}/${DECLINE_RATE}/selection_${S}/${AF_TAG}/model_${MSPRIME_MODEL}"

    OUTDIR="${RUNBASE}/replicate_${REP}"

    mkdir -p "${OUTDIR}"

    cd "${OUTDIR}"


    echo "----------------------------------------"
    echo "REP=${REP}"
    echo "SEED=${SEED}"
    echo "PWD=${PWD}"
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


    # ========================================================
    # Verify outputs exist
    # ========================================================

    if [[ ! -s "${TREES}" ||
          ! -s "${T1_IDS}" ||
          ! -s "${T2_IDS}" ||
          ! -s "${T3_IDS}" ]]; then

        echo "ERROR: Missing expected SLiM outputs in ${OUTDIR}" >&2

        echo "Expected:"
        echo "  ${TREES}"
        echo "  ${T1_IDS}"
        echo "  ${T2_IDS}"
        echo "  ${T3_IDS}"

        echo
        echo "Files present:"

        ls -lh

        exit 2
    fi


    # ========================================================
    # Recapitate + overlay mutations + calculate statistics
    #
    # Calculates statistics for:
    #
    #   T1
    #   T2
    #   T3
    #
    # including:
    #
    #   pi_t1
    #   pi_t2
    #   pi_t3
    #
    #   delta-pi T1 -> T2
    #   delta-pi T2 -> T3
    #   delta-pi T1 -> T3
    #
    # and corresponding delta-AF statistics.
    # ========================================================

    echo "[PY] recap+stats starting..."


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

done


# ============================================================
# Remove temporary replicate directories
#
# Summary statistics have already been appended to STATS_TSV.
# ============================================================

echo "Removing temporary replicate directories..."

rm -rf "${RUNBASE}"


# ============================================================
# Finished
# ============================================================

echo "============================================================"
echo "ALL DONE"
echo "JOB=${SLURM_JOB_ID}"
echo "NREPS=${NREPS}"
echo "SUMMARY=${STATS_TSV}"
echo "============================================================"