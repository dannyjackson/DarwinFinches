#!/bin/bash
#SBATCH --job-name=gsea_power_emp
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G
#
# One condition per job: 100 iterations x 3 species = 300 fgsea runs.
# RESUMABLE: rerun with the same flags and finished iterations are skipped.
#
# Usage:
#   sbatch run_gsea_power_empirical.sh -N 100k -o 5 -r 1.00 -s 0.1 -f 0.05 [-n 100] [-g GO:0007596]
#   (-s 0.0 runs the null control; -f is then ignored)

set -euo pipefail

N_ITER=100
TERM="GO:0007596"
F0=0.5

while getopts "N:o:r:s:f:n:g:" option; do
    case "${option}" in
        N) NE=${OPTARG} ;;          # directory spelling: 50k, 100k, 200k
        o) GEN=${OPTARG} ;;
        r) DECLINE=${OPTARG} ;;
        s) S=${OPTARG} ;;
        f) F0=${OPTARG} ;;
        n) N_ITER=${OPTARG} ;;
        g) TERM=${OPTARG} ;;
        *) echo "Invalid option" >&2; exit 1 ;;
    esac
done

: "${NE:?-N required}" "${GEN:?-o required}" "${DECLINE:?-r required}" "${S:?-s required}"

PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026

cd "${PROJDIR}"
module load R

Rscript "${SCRIPTDIR}/gsea_power_empirical.R" \
    Ne="${NE}" Gen="${GEN}" decline="${DECLINE}" s="${S}" f0="${F0}" \
    n_iter="${N_ITER}" term="${TERM}" \
    nproc="${SLURM_CPUS_PER_TASK:-1}" \
    out_dir="${PROJDIR}/gsea_empirical_power"
