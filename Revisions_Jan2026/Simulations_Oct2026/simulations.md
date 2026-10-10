# GSEA Simulations

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

# test script
F0=0.5
S=0.0
L=10000
NE_VALUES=(50000)
DECLINES=(1.00)
GEN_TIMES=(10)

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

```
## Model positive selection
```
F0S=(0.05 0.10 0.25 0.5)
NES=(50000 100000 200000)
SELS=(0.01 0.05 0.10 0.2 0.25)
DECLINES=(1.00 0.80)
GEN_TIMES=(5 10)
LENS=(10000)

F0S=(0.05)
NES=(50000)
SELS=(0.25)
DECLINES=(1.00)
GEN_TIMES=(10)
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
```


Rscript ${SCRIPTDIR}/plot.power_analysis.gsea.R