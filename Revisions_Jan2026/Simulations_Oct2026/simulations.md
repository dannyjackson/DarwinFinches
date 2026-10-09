# GSEA Simulations

## Model neutral evolution
```
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/
cd $PROJDIR
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026/

F0=0.5
S=0.0
L=10000
MODEL="hudson"
NE_VALUES=(100000 200000)
DECLINES=(1.00 0.8)
GEN_TIMES=(5)



for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      echo "Submitting: Ne=${NE} decline=${D} gen_time=${G}"

      sbatch ${SCRIPTDIR}/run_slim.power_analysis.sh \
        -s "${S}" \
        -r "${D}" \
        -f "${F0}" \
        -m "${MODEL}" \
        -N "${NE}" \
        -o "${G}" \
        -L "${L}" \
        -n 10000

    done
  done
done

```
## Model positive selection
```
MODEL="hudson"
F0S=(0.05 0.10 0.25 0.5)
# NES=(50000 100000 200000)
SELS=(0.01 0.05 0.10 0.2 0.25)
DECLINES=(1.00 0.80)
GEN_TIMES=(5 10)
LENS=(10000)



for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G in "${GEN_TIMES[@]}"; do
          for L in "${LENS[@]}"; do

            echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} gen=${G} L=${L}"

            sbatch ${SCRIPTDIR}/run_slim.power_analysis.short.sh \
              -s "${S}" \
              -r "${D}" \
              -f "${F0}" \
              -m "${MODEL}" \
              -N "${NE}" \
              -o "${G}" \
              -L "${L}" \
              -n 1000

            # Optional throttle to avoid scheduler spam
            sleep 0.1

          done
        done
      done
    done
  done
done 