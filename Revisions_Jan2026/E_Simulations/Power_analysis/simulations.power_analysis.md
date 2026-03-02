# Modeling sensitivity of our sampling scheme and analytical approach

cd /xdisk/mcnew/finches/dannyjackson/simulations/power_analysis/

mkdir neutral

How much data will this generate?
2.3G per run of 100
5 declines, 4 Ne, 3 gen times.
2.3*1000*5*4*3

Check on progress:
ls */*/*0.0.0.*

# This is intractible as a loop. Rather, run it over a handful of parameters
## NE 10k
### DECLINES 1.00
#### GEN_TIMES 5
```
F0=(0.5)
NE=(10000)
S=(0.0)
D=(1.00)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(10000)
S=(0.0)
D=(1.00)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
### DECLINES 0.80
#### GEN_TIMES 5
```
F0=(0.5)
NE=(10000)
S=(0.0)
D=(0.80)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(10000)
S=(0.0)
D=(0.80)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
## NE 50k
### DECLINES 1.00
#### GEN_TIMES 5
```
F0=(0.5)
NE=(50000)
S=(0.0)
D=(1.00)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(50000)
S=(0.0)
D=(1.00)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
### DECLINES 0.80
#### GEN_TIMES 5
```
F0=(0.5)
NE=(50000)
S=(0.0)
D=(0.80)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(50000)
S=(0.0)
D=(0.80)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
## NE 100k
### DECLINES 1.00
#### GEN_TIMES 5
```
F0=(0.5)
NE=(100000)
S=(0.0)
D=(1.00)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(100000)
S=(0.0)
D=(1.00)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
### DECLINES 0.80
#### GEN_TIMES 5
```
F0=(0.5)
NE=(100000)
S=(0.0)
D=(0.80)
G=(5)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```
#### GEN_TIMES 15
```
F0=(0.5)
NE=(100000)
S=(0.0)
D=(0.80)
G=(15)
L=(10000)
MODEL="hudson"

sbatch run_slim.power_analysis.sh \
    -s "${S}" \
    -r "${D}" \
    -f "${F0}" \
    -m "${MODEL}" \
    -N "${NE}" \
    -o "${G}" \
    -L "${L}" \
    -n 10000
```




F0S=(0.5)
NES=(10000 100000)
SELS=(0.0)
DECLINES=(1.00 0.80)
GEN_TIMES=(5 15)
LENS=(10000)

```
# remove all prior runs of selection

wc -l */*/*0.00.summary_stats.tsv
wc -l */*/*0.01.summary_stats.tsv
wc -l */*/*0.05.summary_stats.tsv
wc -l */*/*0.10.summary_stats.tsv
wc -l */*/*0.20.summary_stats.tsv


rm */*/*0.00.summary_stats.tsv
rm */*/*0.01.summary_stats.tsv
rm */*/*0.05.summary_stats.tsv
rm */*/*0.10.summary_stats.tsv
rm */*/*0.20.summary_stats.tsv
rm */*/*0.30.summary_stats.tsv
rm */*/*0.40.summary_stats.tsv
rm */*/*0.50.summary_stats.tsv



# Reduced parameter space
F0S=(0.5)
NES=(10000 100000)
SELS=(0.00 0.01 0.05 0.10 0.20 0.30 0.40 0.50)
DECLINES=(1.00 0.80)
GEN_TIMES=(5 15)
LENS=(10000)

for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G in "${GEN_TIMES[@]}"; do
          for L in "${LENS[@]}"; do

            echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} gen=${G} L=${L}"

            sbatch run_slim.power_analysis.sh \
              -s "${S}" \
              -r "${D}" \
              -f "${F0}" \
              -m "${MODEL}" \
              -N "${NE}" \
              -o "${G}" \
              -L "${L}" \
              -n 102

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

Rscript plot.power_analysis.r



10000	5	0.8	0	0.01
10000	5	0.8	0	0.05
10000	5	0.8	0	0.1
10000	5	0.8	0	0.2
10000	5	0.8	0	0.3
10000	5	0.8	0	0.4
10000	5	0.8	0	0.5

