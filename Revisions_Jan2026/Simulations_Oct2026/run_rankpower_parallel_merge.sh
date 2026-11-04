#!/bin/bash
#SBATCH --job-name=gsea_par_merge
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.err
#SBATCH --time=00:10:00
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#
# Concatenates output/parallel3_parts/*.tsv into
# output/fst_rank_by_selected_windows.parallel3.nnull10000.tsv and lists any
# condition without a part file (rerun those tasks before trusting the merge).

set -euo pipefail
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
cd "$PROJDIR"
module load R
Rscript "${SCRIPTDIR}/D1_RankPower.top5.parallel.R" merge
