#!/bin/bash
#SBATCH --job-name=slim_sweep
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%A_%a.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#
# Measured peak memory per replicate: ~0.7 GB at Ne=200k (SLiM step).
# Measured time per replicate on local disk: ~12 s (Ne=50k), ~30 s (100k),
# ~50 s (200k), plus ~40 s per extra msprime attempt at 200k. With -n 100,
# a task takes ~0.5 h (50k) to ~3 h (200k), well inside 12 h.
#
# RESUMABLE: replicates already in the summary TSV are skipped. If tasks time
# out or fail, resubmit the SAME --array indices with the SAME flags; they
# continue where they stopped and never duplicate a replicate.
#
# One parameter combination per submission. Replicates per task = -n.
# Submit as an array to spread replicates over many tasks, e.g.
#   sbatch --array=1-100 run_slim.power_analysis.sh ... -n 100   # 10,000 replicates
# Task k runs replicates (k-1)*n+1 .. k*n, so replicate IDs never overlap.
#
# Per replicate:
#   1. make_initial_state.py      msprime genealogy of all Ne diploids + neutral
#                                 mutations; picks an existing variant at ~f0
#   2. SLiM                       loads that population, samples T1, applies
#                                 selection (s may be 0) and decline, samples T2
#   3. compute_stats.power_analysis.py   FST, pi, dxy, Tajima's D, selected-site freqs
#
# Rough cost at Ne=100k, L=10 kb (single core): ~10-25 s for step 1
# (more when a retry is needed), ~5 s SLiM, ~3 s stats. Roughly double at Ne=200k.

set -euo pipefail

NREPS=100
TOL=""

while getopts "s:r:f:m:N:o:L:n:t:" option; do
    case "${option}" in
        s) S=${OPTARG} ;;
        r) DECLINE_RATE=${OPTARG} ;;
        f) F0=${OPTARG} ;;
        m) MSPRIME_MODEL=${OPTARG} ;;   # kept for compatibility; ancestry is DTWF->Hudson
        N) NE=${OPTARG} ;;
        o) OFFSET=${OPTARG} ;;
        L) SEG_LEN=${OPTARG} ;;
        n) NREPS=${OPTARG} ;;
        t) TOL=${OPTARG} ;;             # optional |freq - f0| tolerance for the chosen variant
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1 ;;
    esac
done

: "${S:?ERROR: -s <sel_s> is required}"
: "${DECLINE_RATE:?ERROR: -r <decline_rate> is required}"
F0="${F0:-0.10}"
MSPRIME_MODEL="${MSPRIME_MODEL:-hudson}"
NE="${NE:-10000}"
OFFSET="${OFFSET:-10}"
SEG_LEN="${SEG_LEN:-10000}"

TASK="${SLURM_ARRAY_TASK_ID:-1}"
REP_START=$(( (TASK - 1) * NREPS + 1 ))
REP_END=$(( TASK * NREPS ))

# ----------------------------
# Paths
# ----------------------------
SIMDIR="/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis"
SCRIPTDIR="${SIMDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026/"
SLIM_BIN="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5/bin/slim"
SLIM_SCRIPT="${SCRIPTDIR}/simulation_replicates.power_analysis.slim"
PREP_PY="${SCRIPTDIR}/make_initial_state.py"
STATS_PY="${SCRIPTDIR}/compute_stats.power_analysis.py"
PYTHON="/xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py/bin/python"

NE_TAG="Ne_${NE/000/}k"
GEN_TAG="Gen_${OFFSET}"
AF_TAG="AF_${F0//./_}"

OUTBASE="${SIMDIR}/${NE_TAG}/${GEN_TAG}"
mkdir -p "${OUTBASE}/logs"

# Same file name as before, so the plotting / RankPower scripts still find it
STATS_TSV="${OUTBASE}/d${DECLINE_RATE}_s${S}_af${F0}.summary_stats.tsv"
LOCKFILE="${STATS_TSV}.lock"

# Temporary replicate dirs are unique to this parameter combination AND task,
# so concurrent jobs (e.g. different f0) can never overwrite or delete each other's files.
# Work on the compute node's local disk, not /xdisk: each replicate writes and
# reads two ~15-30 MB tree files, and 100 tasks doing that on the shared
# filesystem at once is the main slowdown. Only the one-line result row and the
# trajectories are written to /xdisk.
LOCAL_SCRATCH="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}"
RUNBASE="${LOCAL_SCRATCH}/${USER}/slim_d${DECLINE_RATE}_s${S}_${AF_TAG}_${NE_TAG}_${GEN_TAG}_task${TASK}_job${SLURM_JOB_ID:-local}"
trap 'rm -rf "${RUNBASE}"' EXIT   # clean node-local disk even on timeout/failure

# One thread per tool, matching --cpus-per-task=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# Trajectories go straight to /xdisk (tiny), so they survive a timeout
mkdir -p "${OUTBASE}/trajectories"
TRAJ_FILE="${OUTBASE}/trajectories/d${DECLINE_RATE}_s${S}_af${F0}_task${TASK}.txt"

# Seed base unique to the parameter combination (replicate ID is added below).
# Kept < 2^32 including make_initial_state.py's per-attempt offsets.
PARAM_KEY="Ne${NE}_G${OFFSET}_d${DECLINE_RATE}_s${S}_f${F0}_L${SEG_LEN}"
SEED_BASE=$(( ( $(printf '%s' "${PARAM_KEY}" | cksum | cut -d' ' -f1) % 40000 ) * 100000 ))
if (( REP_END >= 100000 )); then
    echo "ERROR: replicate IDs must stay below 100000 for unique seeds" >&2; exit 1
fi

echo "JOB=${SLURM_JOB_ID:-local} TASK=${TASK} REPS=${REP_START}-${REP_END}"
echo "S=${S} DECLINE_RATE=${DECLINE_RATE} F0=${F0} NE=${NE} OFFSET=${OFFSET} L=${SEG_LEN}"
echo "PARAM_KEY=${PARAM_KEY} SEED_BASE=${SEED_BASE}"
echo "STATS_TSV=${STATS_TSV}"

# Replicates already finished for this condition (column 1 of the TSV)
declare -A DONE=()
if [[ -s "${STATS_TSV}" ]]; then
    while read -r r; do DONE["$r"]=1; done < <(awk -F'\t' 'NR>1{print $1}' "${STATS_TSV}")
fi
N_SKIP=0
for REP in $(seq "${REP_START}" "${REP_END}"); do [[ -n "${DONE[$REP]:-}" ]] && N_SKIP=$((N_SKIP+1)); done
echo "Already done in this task's range: ${N_SKIP} of $((REP_END-REP_START+1))"

TOL_ARG=()
[[ -n "${TOL}" ]] && TOL_ARG=(--tol "${TOL}")

for REP in $(seq "${REP_START}" "${REP_END}"); do
    if [[ -n "${DONE[$REP]:-}" ]]; then continue; fi
    SEED=$(( SEED_BASE + REP ))
    OUTDIR="${RUNBASE}/replicate_${REP}"
    mkdir -p "${OUTDIR}"
    cd "${OUTDIR}"

    echo "---- REP=${REP} SEED=${SEED} start=$(date +%T)"
    t0=${SECONDS}

    # 1. Ancestry-aware starting population
    "${PYTHON}" "${PREP_PY}" \
        --Ne "${NE}" --L "${SEG_LEN}" --f0 "${F0}" ${TOL_ARG[@]+"${TOL_ARG[@]}"} \
        --mu 2.04e-9 --recomb 1e-8 --seed "${SEED}" \
        --out_trees "init_run${REP}.trees" \
        --out_site  "site_run${REP}.txt" \
        > "init_run${REP}.log" 2>&1
    t1=${SECONDS}

    # 2. Forward simulation
    "${SLIM_BIN}" \
        -d run_id="${REP}" \
        -d seed="${SEED}" \
        -d sel_s="${S}" \
        -d decline_rate="${DECLINE_RATE}" \
        -d f0="${F0}" \
        -d Ne="${NE}" \
        -d L="${SEG_LEN}" \
        -d T2_OFFSET="${OFFSET}" \
        -d "init_trees='init_run${REP}.trees'" \
        -d "site_file='site_run${REP}.txt'" \
        "${SLIM_SCRIPT}" > "run${REP}.log" 2>&1
    t2=${SECONDS}

    for f in "simulation_run${REP}.trees" "sample_t1_run${REP}.ids.txt" "sample_t2_run${REP}.ids.txt"; do
        [[ -s "$f" ]] || { echo "ERROR: missing $f in ${OUTDIR}"; tail -20 "run${REP}.log"; exit 2; }
    done

    # 3. Statistics -> one-row TSV, then append under a lock
    "${PYTHON}" "${STATS_PY}" \
        --trees "simulation_run${REP}.trees" \
        --t1_ids "sample_t1_run${REP}.ids.txt" \
        --t2_ids "sample_t2_run${REP}.ids.txt" \
        --site_file "site_run${REP}.txt" \
        --Ne "${NE}" --mu 2.04e-9 --recomb 1e-8 \
        --offset "${OFFSET}" --seed "${SEED}" \
        --rep "${REP}" --sel_s "${S}" --decline_rate "${DECLINE_RATE}" --f0 "${F0}" \
        --tsv_out "row.tsv" \
        > "stats_run${REP}.log" 2>&1
    t3=${SECONDS}

    (
        flock -x 9
        if [[ ! -s "${STATS_TSV}" ]]; then
            cat row.tsv >> "${STATS_TSV}"
        else
            tail -n +2 row.tsv >> "${STATS_TSV}"
        fi
    ) 9> "${LOCKFILE}"
    t4=${SECONDS}
    echo "     timing_s init=$((t1-t0)) slim=$((t2-t1)) stats=$((t3-t2)) append=$((t4-t3)) $(grep -c '^attempt' "init_run${REP}.log") msprime attempt(s)"

    # Keep the selected-allele trajectory, drop the large files
    grep "selected" "run${REP}.log" | sed "s/^/rep${REP}\t/" >> "${TRAJ_FILE}" || true
    cd "${RUNBASE}"
    rm -rf "${OUTDIR}"
done

rm -rf "${RUNBASE}"

echo "ALL DONE (JOB=${SLURM_JOB_ID:-local}, TASK=${TASK})"
