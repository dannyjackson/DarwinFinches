#!/usr/bin/env bash
# Top up the positive-selection grid: bring every condition to replicates
# 1..TARGET_REPS, submitting only what is missing and only as many array tasks
# as your submit limit allows. Run it again whenever jobs finish; it carries on
# where it left off and says when everything is done.
#
# Run on the login node with bash (not sbatch):
#   DRYRUN=1 bash topup_selected.sh    # preview
#   bash topup_selected.sh             # submit
#
# How it stays safe while earlier jobs are still running:
#   - It reads the full sbatch command of every queued/running array job
#     (sacct SubmitLine) and SKIPS any condition that still has jobs in flight,
#     so the same replicate can never run twice at once.
#   - Replicate IDs and seeds depend only on the condition and the replicate
#     number, and run_slim.power_analysis.sh skips replicates already in the
#     summary file. So a condition can be finished with a different -n than it
#     was started with: the replicates (and results) are identical.
#
# Fewer, larger tasks: replicates per task scale with Ne so each task takes a
# few hours (measured on local disk: ~12 s/rep at 50k, ~30 s at 100k,
# ~50-90 s at 200k). That cuts the grid from 3,600 tasks to ~840.

set -euo pipefail

PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
RUN_SCRIPT=${SCRIPTDIR}/run_slim.power_analysis.sh
cd "${PROJDIR}"

TARGET_REPS=1000          # replicates 1..1000 per selected condition

# Replicates per array task, by Ne (keep each task well under the 12 h limit)
declare -A REPS_PER_TASK=( [50000]=1000 [100000]=500 [200000]=250 )
DEFAULT_REPS_PER_TASK=100

# Max queued + running array tasks you are allowed. Auto-detected from your
# QOS/association if possible; otherwise this default. Override: QOS_LIMIT=500
QOS_LIMIT_DEFAULT=1000
SAFETY=10                 # leave a little headroom under the limit

NES=(50000 100000 200000)
SELS=(0.01 0.02 0.05 0.1 0.2 0.5)
F0S=(0.05 0.1 0.2 0.4 0.6)
DECLINES=(1.00 0.8)        # same spelling as the original submissions
GEN_TIMES=(5 10)
L=10000

# ------------------------------------------------------------
# Submit limit and current usage
# ------------------------------------------------------------
detect_limit() {
  local a q
  a=$(sacctmgr -n -P show assoc user="${USER}" format=maxsubmit 2>/dev/null | awk -F'|' '$1+0>0{print $1+0}' | sort -n | head -1)
  q=$(sacctmgr -n -P show qos format=maxsubmitpu 2>/dev/null | awk -F'|' '$1+0>0{print $1+0}' | sort -n | head -1)
  printf '%s\n' "${a}" "${q}" | awk 'NF && $1>0' | sort -n | head -1
}
if [[ -n "${QOS_LIMIT:-}" ]]; then
  LIMIT=${QOS_LIMIT}
else
  LIMIT=$(detect_limit || true)
  LIMIT=${LIMIT:-${QOS_LIMIT_DEFAULT}}
fi
ACTIVE_TASKS=$(squeue -u "${USER}" -h -r -t PD,R,CF,CG 2>/dev/null | wc -l)
BUDGET=$(( LIMIT - ACTIVE_TASKS - SAFETY ))
echo "Submit limit ${LIMIT}; ${ACTIVE_TASKS} task(s) queued/running; room for ${BUDGET} more"

# ------------------------------------------------------------
# Conditions that still have jobs in flight (from the sbatch command line)
# Key: Ne|gen|decline|s|f0 with numbers normalised (0.80 == 0.8, 1.00 == 1)
# ------------------------------------------------------------
norm() { awk -v x="$1" 'BEGIN{printf "%g", x}'; }
key()  { echo "$(norm "$1")|$(norm "$2")|$(norm "$3")|$(norm "$4")|$(norm "$5")"; }

declare -A BUSY=()
UNKNOWN_ACTIVE=0
ARRAY_IDS=$(squeue -u "${USER}" -h -t PD,R,CF,CG -o "%F" 2>/dev/null | sort -u)
for J in ${ARRAY_IDS}; do
  LINE=$(sacct -j "${J}" -X -n -P --format=SubmitLine 2>/dev/null | head -1 || true)
  if [[ -z "${LINE}" || "${LINE}" != *run_slim.power_analysis.sh* ]]; then
    # Not a simulation job, or SubmitLine unavailable
    [[ "$(squeue -j "${J}" -h -o %j 2>/dev/null | head -1)" == slim_sweep ]] && UNKNOWN_ACTIVE=1
    continue
  fi
  read -r s r f n o < <(awk '{
      for (i=1;i<NF;i++) {
        if ($i=="-s") s=$(i+1); if ($i=="-r") r=$(i+1); if ($i=="-f") f=$(i+1)
        if ($i=="-N") n=$(i+1); if ($i=="-o") o=$(i+1)
      }
      gsub(/"/,"",s); gsub(/"/,"",r); gsub(/"/,"",f); gsub(/"/,"",n); gsub(/"/,"",o)
      print s, r, f, n, o }' <<< "${LINE}")
  BUSY["$(key "${n}" "${o}" "${r}" "${s}" "${f}")"]=1
done
if (( UNKNOWN_ACTIVE )) && [[ "${FORCE:-0}" != 1 ]]; then
  echo "ERROR: some slim_sweep jobs are running but their command lines can't be read"
  echo "(sacct SubmitLine unavailable), so in-flight conditions can't be identified."
  echo "Wait for them to finish, or rerun with FORCE=1 if you're sure there's no overlap."
  exit 1
fi
echo "Conditions with jobs in flight (skipped this round): ${#BUSY[@]}"

to_ranges() {
  awk 'BEGIN{ORS=""}
       { if (NR==1) {s=$1; p=$1; next}
         if ($1==p+1) {p=$1; next}
         print (s==p ? s : s"-"p) ","; s=$1; p=$1 }
       END{ if (NR>0) print (s==p ? s : s"-"p) }'
}

N_DONE=0; N_BUSY=0; N_SUBMITTED=0; N_WAITING=0; TASKS_SUBMITTED=0; STOP=0

for F0 in "${F0S[@]}"; do
 for NE in "${NES[@]}"; do
  for S in "${SELS[@]}"; do
   for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      LABEL="f0=${F0} Ne=${NE} s=${S} d=${D} gen=${G}"
      NE_TAG="Ne_${NE/000/}k"
      TSV="${PROJDIR}/${NE_TAG}/Gen_${G}/d${D}_s${S}_af${F0}.summary_stats.tsv"

      if [[ -n "${BUSY[$(key "${NE}" "${G}" "${D}" "${S}" "${F0}")]:-}" ]]; then
        N_BUSY=$(( N_BUSY + 1 )); continue
      fi

      # Which of replicates 1..TARGET_REPS are already done
      if [[ -s "${TSV}" ]]; then
        DONE_LIST=$(awk -F'\t' -v t="${TARGET_REPS}" 'NR>1 && $1+0>=1 && $1+0<=t {print $1+0}' "${TSV}" | sort -n -u)
      else
        DONE_LIST=""
      fi
      N_HAVE=$(grep -c . <<< "${DONE_LIST}" || true)
      if (( N_HAVE >= TARGET_REPS )); then
        N_DONE=$(( N_DONE + 1 )); continue
      fi

      if (( STOP )); then N_WAITING=$(( N_WAITING + 1 )); continue; fi

      NPT=${REPS_PER_TASK[$NE]:-${DEFAULT_REPS_PER_TASK}}
      NTASK_MAX=$(( (TARGET_REPS + NPT - 1) / NPT ))

      # Task indices (for this -n) that still have missing replicates
      PICKED=$(awk -v n="${NPT}" -v m="${NTASK_MAX}" '
        { done[int(($1-1)/n)+1]++ }
        END { for (t=1; t<=m; t++) {
                last = (t*n > '"${TARGET_REPS}"') ? '"${TARGET_REPS}"' : t*n
                need = last - (t-1)*n
                if (done[t]+0 < need) print t } }' <<< "${DONE_LIST}")
      NT=$(grep -c . <<< "${PICKED}" || true)
      (( NT == 0 )) && { N_DONE=$(( N_DONE + 1 )); continue; }

      if (( NT > BUDGET )); then
        STOP=1; N_WAITING=$(( N_WAITING + 1 )); continue
      fi

      ARRAY=$(to_ranges <<< "${PICKED}")
      echo "${LABEL}: ${N_HAVE}/${TARGET_REPS} done -> --array=${ARRAY} -n ${NPT}"
      CMD=(sbatch --array="${ARRAY}" "${RUN_SCRIPT}"
           -s "${S}" -r "${D}" -f "${F0}" -N "${NE}" -o "${G}" -L "${L}" -n "${NPT}")

      if [[ "${DRYRUN:-0}" == 1 ]]; then
        echo "    DRYRUN: ${CMD[*]}"
      else
        if ! "${CMD[@]}"; then
          echo "    sbatch refused (probably the submit limit); stopping this round"
          STOP=1; N_WAITING=$(( N_WAITING + 1 )); continue
        fi
        sleep 0.1
      fi
      BUDGET=$(( BUDGET - NT ))
      TASKS_SUBMITTED=$(( TASKS_SUBMITTED + NT ))
      N_SUBMITTED=$(( N_SUBMITTED + 1 ))

    done
   done
  done
 done
done

TOTAL=$(( ${#F0S[@]} * ${#NES[@]} * ${#SELS[@]} * ${#DECLINES[@]} * ${#GEN_TIMES[@]} ))
echo "------------------------------------------------------------"
VERB="submitted now"; [[ "${DRYRUN:-0}" == 1 ]] && VERB="would be submitted"
echo "Conditions: ${TOTAL} total | ${N_DONE} complete | ${N_BUSY} with jobs in flight | ${N_SUBMITTED} ${VERB} (${TASKS_SUBMITTED} tasks) | ${N_WAITING} waiting for room"
if (( N_DONE == TOTAL )); then
  echo "ALL CONDITIONS COMPLETE"
elif (( N_WAITING > 0 || N_BUSY > 0 )); then
  echo "Run this script again after some jobs finish."
fi
