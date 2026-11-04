#!/bin/bash
#SBATCH --job-name=gsea_par_arr
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%A_%a.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%A_%a.err
#SBATCH --time=03:00:00
#SBATCH --array=0-47
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# Three-species gene-set rank test (approach 2b, synthetic) as a Slurm array.
# 480 conditions split over 48 tasks = 10 conditions per task; at the observed
# ~3.6 min per condition that is ~40 min per task, so 3 h is a wide margin.
# Each condition writes its own part file; rerunning the array skips finished
# conditions, so a failed task can be resubmitted with --array=<id>.
#
# Submit:  sbatch run_rankpower_parallel_array.sh
# Merge:   sbatch --dependency=afterok:<array job id> run_rankpower_parallel_merge.sh
#          (or run the merge line by hand once squeue shows no tasks left)

set -euo pipefail
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
cd "$PROJDIR"
module load R

N_TASKS=$(( SLURM_ARRAY_TASK_MAX - SLURM_ARRAY_TASK_MIN + 1 ))
Rscript "${SCRIPTDIR}/D1_RankPower.top5.parallel.R" "$SLURM_ARRAY_TASK_ID" "$N_TASKS"
