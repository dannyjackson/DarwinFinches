#!/usr/bin/env bash
# Top up neutral replicate files to a target number of usable rows.
#
# For each condition it:
#   - counts rows with a finite FST (D1_RankPower drops NaN/Inf FST rows)
#   - finds the highest replicate ID already in the file
#   - submits just enough new array tasks, numbered so their replicate IDs
#     (and therefore seeds) start above every ID already used
#
# Replicate IDs for task k are (k-1)*NREPS+1 .. k*NREPS, so new tasks never
# repeat an existing replicate as long as NREPS stays the same.
#
# Run from the project directory. Set DRYRUN=1 to print the sbatch
# commands without submitting.

set -euo pipefail

PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
cd "${PROJDIR}"

TARGET=10000      # usable rows wanted per condition
MARGIN=300        # extra replicates to cover NaN-FST rows and failed tasks
NREPS=100         # replicates per array task (keep equal to the original run)

F0=0.5
S=0.0
L=10000
NE_VALUES=(50000)
DECLINES=(1.00)
GEN_TIMES=(10)

for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      NE_TAG="Ne_${NE/000/}k"
      TSV="${PROJDIR}/${NE_TAG}/Gen_${G}/d${D}_s${S}_af${F0}.summary_stats.tsv"

      if [[ -s "${TSV}" ]]; then
        # column 1 = rep, column 12 = fst_hudson
        read -r USABLE MAXREP < <(awk -F'\t' '
          NR > 1 {
            if ($12 != "nan" && $12 != "inf" && $12 != "-inf" && $12 != "") n++
            if ($1 + 0 > m) m = $1 + 0
          }
          END { print n + 0, m + 0 }' "${TSV}")
      else
        USABLE=0; MAXREP=0
      fi

      NEED=$(( TARGET - USABLE ))
      if (( NEED <= 0 )); then
        echo "Ne=${NE} d=${D} gen=${G}: ${USABLE} usable rows, nothing to do"
        continue
      fi

      NTASKS=$(( (NEED + MARGIN + NREPS - 1) / NREPS ))
      # First task whose replicate IDs start above MAXREP (and above the
      # original 1..100 tasks, which may still be partly running)
      FIRST=$(( (MAXREP + NREPS - 1) / NREPS + 1 ))
      (( FIRST < 101 )) && FIRST=101
      LAST=$(( FIRST + NTASKS - 1 ))

      echo "Ne=${NE} d=${D} gen=${G}: ${USABLE} usable rows (max rep ${MAXREP})," \
           "submitting tasks ${FIRST}-${LAST} = reps $(( (FIRST-1)*NREPS+1 ))-$(( LAST*NREPS ))"

      CMD=(sbatch --array="${FIRST}-${LAST}"
           "${SCRIPTDIR}/run_slim.power_analysis.sh"
           -s "${S}" -r "${D}" -f "${F0}" -N "${NE}" -o "${G}" -L "${L}" -n "${NREPS}")

      if [[ "${DRYRUN:-0}" == 1 ]]; then
        echo "  DRYRUN: ${CMD[*]}"
      else
        "${CMD[@]}"
      fi

    done
  done
done
