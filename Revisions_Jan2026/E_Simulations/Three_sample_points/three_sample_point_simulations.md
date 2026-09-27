# Modeling sensitivity of possible extended sampling schemes and analytical approaches

cd /xdisk/mcnew/dannyjackson/

micromamba create -n slim5

micromamba install -p /xdisk/mcnew/dannyjackson/.local/share/mamba/envs/slim5 \
    -c conda-forge -c bioconda \
    --force-reinstall slim

micromamba create \
    -p /xdisk/mcnew/dannyjackson/.local/share/mamba/envs/recap_py \
    -c conda-forge \
    python=3.12 \
    numpy \
    tskit \
    msprime \
    pyslim

cd /xdisk/mcnew/finches/dannyjackson/simulations/three_sample_points

# revise sh script to direct it to simulation_replicates.three_sample_points.slim
# revise to also compute delta AF instead of fst
# compute delta pi as well




```

MODEL="hudson"

# F0S=(0.05 0.25 0.50 0.75 0.95)
# NES=(10000 50000 100000 150000)
# SELS=(0.00 0.01 0.05 0.10 0.20 0.5)
# DECLINES=(1.00 0.95 0.90 0.85 0.80)

# Generations between t1 and t2
# GEN_TIMES_T1_T2=(5 10 15)

# Generations between t2 and t3
# GEN_TIMES_T2_T3=(5 10 15)

LENS=(10000)

# Test values

F0S=(0.05)
NES=(150000)
SELS=(0.5)
DECLINES=(1.00)

# Generations between t1 and t2
GEN_TIMES_T1_T2=(15)

# Generations between t2 and t3
GEN_TIMES_T2_T3=(5)

for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G12 in "${GEN_TIMES_T1_T2[@]}"; do
          for G23 in "${GEN_TIMES_T2_T3[@]}"; do
            for L in "${LENS[@]}"; do

              # SLiM offsets are both measured from t1
              T2_OFFSET="${G12}"
              T3_OFFSET=$((G12 + G23))

              echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} t1-t2=${G12} t2-t3=${G23} T2_OFFSET=${T2_OFFSET} T3_OFFSET=${T3_OFFSET} L=${L}"

              sbatch run_slim.three_sample_points.sh \
                -s "${S}" \
                -r "${D}" \
                -f "${F0}" \
                -m "${MODEL}" \
                -N "${NE}" \
                -o "${T2_OFFSET}" \
                -O "${T3_OFFSET}" \
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
done

```
