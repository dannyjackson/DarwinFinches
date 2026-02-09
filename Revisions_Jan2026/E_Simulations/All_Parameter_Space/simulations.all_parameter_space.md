# Modeling sensitivity of our sampling scheme and analytical approach

cd /xdisk/mcnew/finches/dannyjackson/simulations/all_parameter_space

```

MODEL="hudson"

# DECLINES=(0.00 0.05 0.10 0.15 0.20)



F0S=(0.05 0.25 0.50 0.75 0.95)
NES=(10000 50000 100000 150000)
# SELS=(0.00 0.01 0.05 0.10 0.20)
SELS=(0.5)
DECLINES=(1.00 0.95 0.90 0.85 0.80)
GEN_TIMES=(5 10 15)
LENS=(10000)



for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G in "${GEN_TIMES[@]}"; do
          for L in "${LENS[@]}"; do

            echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} gen=${G} L=${L}"

            sbatch run_slim.all_parameter_space.sh \
              -s "${S}" \
              -r "${D}" \
              -f "${F0}" \
              -m "${MODEL}" \
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


## Plot the results
module load micromamba
micromamba activate r_elgato

Rscript plot.all_parameter_space.r
