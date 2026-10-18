#!/usr/bin/env bash
# Top up every neutral (s = 0.0) condition to at least TARGET usable rows.
#
# Run on the login node with bash (not sbatch):
#   DRYRUN=1 bash topup_neutral.sh     # preview what would be submitted
#   bash topup_neutral.sh              # submit
#
# For each condition it:
#   1. reads the replicate IDs already in the summary file (missing file = 0 rows)
#   2. counts usable rows (finite FST; D1_RankPower drops the rest)
#   3. walks array task indices 1, 2, 3, ... and picks the INCOMPLETE ones until
#      the expected usable total reaches TARGET + MARGIN
#   4. submits only those task indices as one array
#
# Because run_slim.power_analysis.sh skips replicates already in the file,
# re-running an incomplete task only does its missing replicates, with the same
# seeds a complete first run would have used. Gaps left by timed-out tasks are
# filled first, so a condition ends up as close as possible to replicates
# 1..10,000 of an uninterrupted run.
#
# Safety: two copies of the same task running at once would duplicate
# replicates, so the script refuses to submit while any slim_sweep jobs are
# queued or running (override with FORCE=1 only if you are sure they are for
# other conditions).

set -euo pipefail

PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
RUN_SCRIPT=${SCRIPTDIR}/run_slim.power_analysis.sh
cd "${PROJDIR}"

TARGET=10000      # usable rows wanted per condition
MARGIN=100        # small cushion for NaN-FST rows and failures
NREPS=100         # replicates per array task: must match the original runs
# Your QOS rejected submissions at roughly 1,000 queued array tasks
# ("QOSMaxSubmitJobPerUserLimit"). Check the exact limit with:
#   sacctmgr show assoc user=$USER format=user,maxsubmit
#   sacctmgr show qos format=name,maxsubmitpu
# The script stops adding tasks once this many are submitted in one go;
# just run it again after those finish and it continues where it left off.
MAX_TOTAL_TASKS=${MAX_TOTAL_TASKS:-900}

F0=0.5
S=0.0
L=10000
NE_VALUES=(10000 50000 100000 200000)
DECLINES=(1.00 0.8)      # keep the same spelling as the original runs (0.8, not 0.80)
GEN_TIMES=(5 10)

# ------------------------------------------------------------
# Refuse to run while simulation jobs are still queued/running
# ------------------------------------------------------------
if [[ "${FORCE:-0}" != 1 ]]; then
  ACTIVE=$(squeue -u "${USER}" -h -n slim_sweep -r 2>/dev/null | wc -l || echo 0)
  if (( ACTIVE > 0 )); then
    echo "ERROR: ${ACTIVE} slim_sweep array task(s) still queued or running."
    echo "Wait for them to finish (or scancel them) before topping up,"
    echo "otherwise the same replicates could run twice. FORCE=1 overrides."
    exit 1
  fi
fi

# Compress a sorted list of integers into SLURM array syntax: 1-5,8,10-12
to_ranges() {
  awk 'BEGIN{ORS=""}
       { if (NR==1) {s=$1; p=$1; next}
         if ($1==p+1) {p=$1; next}
         print (s==p ? s : s"-"p) ","; s=$1; p=$1 }
       END{ if (NR>0) print (s==p ? s : s"-"p) }'
}

TOTAL_TASKS=0
DEFERRED=0

for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      NE_TAG="Ne_${NE/000/}k"
      TSV="${PROJDIR}/${NE_TAG}/Gen_${G}/d${D}_s${S}_af${F0}.summary_stats.tsv"
      LABEL="Ne=${NE} d=${D} gen=${G}"

      # Unique replicate IDs, usable flag per ID (col 1 = rep, col 12 = fst_hudson)
      if [[ -s "${TSV}" ]]; then
        read -r ROWS UNIQUE USABLE DUPS < <(awk -F'\t' '
          NR > 1 {
            rows++
            if (seen[$1]++) { dups++; next }
            uniq++
            if ($12 != "nan" && $12 != "inf" && $12 != "-inf" && $12 != "") ok++
          }
          END { print rows+0, uniq+0, ok+0, dups+0 }' "${TSV}")
        DONE_LIST=$(awk -F'\t' 'NR>1{print $1}' "${TSV}" | sort -n -u)
      else
        ROWS=0; UNIQUE=0; USABLE=0; DUPS=0; DONE_LIST=""
      fi

      (( DUPS > 0 )) && echo "WARNING ${LABEL}: ${DUPS} duplicated replicate IDs in ${TSV} (counted once)"

      if (( USABLE >= TARGET )); then
        echo "${LABEL}: ${USABLE} usable rows - done"
        continue
      fi

      # Fraction of replicates that give a usable FST (for projecting)
      if (( UNIQUE > 0 )); then
        FRAC=$(awk -v u="${USABLE}" -v n="${UNIQUE}" 'BEGIN{printf "%.6f", u/n}')
      else
        FRAC=1
      fi

      # Replicates done per task index
      declare -A DONE_PER_TASK=()
      if [[ -n "${DONE_LIST}" ]]; then
        while read -r t c; do DONE_PER_TASK[$t]=$c; done < <(
          awk -v n="${NREPS}" '{print int(($1-1)/n)+1}' <<< "${DONE_LIST}" | sort -n | uniq -c | awk '{print $2, $1}')
      fi

      # Walk tasks in order, picking incomplete ones until the projection is enough
      NEEDED=$(( TARGET + MARGIN ))
      PROJECTED=${USABLE}
      PICKED=()
      t=0
      while awk -v p="${PROJECTED}" -v n="${NEEDED}" 'BEGIN{exit !(p < n)}'; do
        t=$(( t + 1 ))
        done_here=${DONE_PER_TASK[$t]:-0}
        missing=$(( NREPS - done_here ))
        (( missing <= 0 )) && continue
        PICKED+=("${t}")
        PROJECTED=$(awk -v p="${PROJECTED}" -v m="${missing}" -v f="${FRAC}" 'BEGIN{printf "%.1f", p + m*f}')
        if (( t > 100000 / NREPS )); then
          echo "ERROR ${LABEL}: would need replicate IDs >= 100000 (seed range); stopping"; exit 1
        fi
      done
      unset DONE_PER_TASK

      NT=${#PICKED[@]}
      ROOM=$(( MAX_TOTAL_TASKS - TOTAL_TASKS ))
      if (( ROOM <= 0 )); then
        echo "${LABEL}: needs ${NT} task(s) - deferred (MAX_TOTAL_TASKS=${MAX_TOTAL_TASKS} reached; run again later)"
        DEFERRED=1
        continue
      fi
      if (( NT > ROOM )); then
        echo "${LABEL}: needs ${NT} task(s), submitting the first ${ROOM} now (task cap); run again later for the rest"
        PICKED=("${PICKED[@]:0:${ROOM}}")
        NT=${ROOM}
        DEFERRED=1
      fi
      ARRAY=$(printf '%s\n' "${PICKED[@]}" | to_ranges)
      TOTAL_TASKS=$(( TOTAL_TASKS + NT ))

      echo "${LABEL}: ${USABLE} usable of ${UNIQUE} reps -> ${NT} task(s), projected ~${PROJECTED%.*} usable"
      echo "    --array=${ARRAY}"

      CMD=(sbatch --array="${ARRAY}" "${RUN_SCRIPT}"
           -s "${S}" -r "${D}" -f "${F0}" -N "${NE}" -o "${G}" -L "${L}" -n "${NREPS}")

      if [[ "${DRYRUN:-0}" == 1 ]]; then
        echo "    DRYRUN: ${CMD[*]}"
      else
        "${CMD[@]}"
      fi

    done
  done
done

if [[ "${DRYRUN:-0}" == 1 ]]; then
  echo "Total array tasks that would be submitted: ${TOTAL_TASKS}"
else
  echo "Total array tasks submitted: ${TOTAL_TASKS}"
fi
if (( DEFERRED )); then
  echo "Some conditions were deferred by the task cap: run this script again once these jobs finish."
fi
