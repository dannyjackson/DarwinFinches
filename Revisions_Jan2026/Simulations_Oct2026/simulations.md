# Simulations

## Set up working environment
```
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026/
cd $SCRIPTDIR 
git pull
cd $PROJDIR
mkdir -p /xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs
```
## Model neutral evolution
```
F0=0.5
S=0.0
L=10000
NE_VALUES=(100000 50000 200000)
DECLINES=(1.00 0.8)
GEN_TIMES=(5 10)

for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      echo "Submitting: Ne=${NE} decline=${D} gen_time=${G}"

      sbatch --array=1-100 \
        ${SCRIPTDIR}/run_slim.power_analysis.sh \
        -s "${S}" \
        -r "${D}" \
        -f "${F0}" \
        -N "${NE}" \
        -o "${G}" \
        -L "${L}" \
         -n 100

    done
  done
done

# topup neutral
DRYRUN=1 bash ${SCRIPTDIR}/topup_neutral.sh   # preview: prints the sbatch line only

bash ${SCRIPTDIR}/topup_neutral.sh

scancel 24232089
squeue -u dannyjackson          # wait until nothing slim_sweep is listed
find Ne_* -name "*.lock" -delete
bash ${SCRIPTDIR}/topup_neutral.sh

```
## Model positive selection
```
NES=(50000 100000 200000)
SELS=(0.01 0.02 0.05 0.1 0.2 0.5)    # 1-2-5 series
SELS=(0.3 0.4)    # 1-2-5 series

F0S=(0.05 0.1 0.2 0.4 0.6)               # doubling
DECLINES=(1.00 0.8)
GEN_TIMES=(5 10)
LENS=(10000)

for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G in "${GEN_TIMES[@]}"; do
          for L in "${LENS[@]}"; do

            echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} gen=${G} L=${L}"

            sbatch --array=1-10 \
               ${SCRIPTDIR}/run_slim.power_analysis.sh \
              -s "${S}" \
              -r "${D}" \
              -f "${F0}" \
              -N "${NE}" \
              -o "${G}" \
              -L "${L}" \
              -n 100

            # Optional throttle to avoid scheduler spam
            sleep 0.1

          done
        done
      done
    done
  done
done 


bash ${SCRIPTDIR}/topup_selected.sh

```

## Model GSEA 

Rscript ${SCRIPTDIR}/plot.power_analysis.gsea.R


```
#!/bin/bash
#SBATCH --job-name=simulation_gsea_plots
#SBATCH --output=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.out
#SBATCH --error=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/output/logs/%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=1
#SBATCH --account=mcnew
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=8G


PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026/
cd $PROJDIR

module load R

Rscript ${SCRIPTDIR}/D1_RankPower.top5.test.R
```
### Modify plots
```
Rscript ${SCRIPTDIR}/plot_power_curves.R
```
## Model outlier analysis
```
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026/
cd $SCRIPTDIR 
git pull
cd $PROJDIR

module load R

Rscript ${SCRIPTDIR}/outlier_power_curves.R
Rscript ${SCRIPTDIR}/plot_false_positive_curves.R
Rscript ${SCRIPTDIR}/plot_ppv_curves.R
```
## Generate multipanel plot of all analyses
```
Rscript ${SCRIPTDIR}/multipanel_power_figure.R

```
## Generate multipanel plot of all analyses using all variations of N
```
cd /xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis

# sweeps 1-100 selected windows -> output/fst_outlier_power_nsel.tsv
Rscript ${SCRIPTDIR}/outlier_power_nsel.R

# In parallel_outlier_overlap.R change line 58 to
# PARALLEL_VALUES <- c(1, 2, 5, 10, 20, 30, 50, 100)
Rscript ${SCRIPTDIR}/parallel_outlier_overlap.R

Rscript ${SCRIPTDIR}/plot_gsea_parallel_nsel.R      # writes summary_nsel.tsv for D

Rscript ${SCRIPTDIR}/multipanel_power_figure_nsel.R
```