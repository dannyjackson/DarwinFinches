# full script


Run each with 1000 replicates

Selection 0.0 0.25 0.5 0.75 1.0
Pop decline 0.05 0.10 0.25 0.5





cd /xdisk/mcnew/finches/dannyjackson/simulations/

wget --output-document /dev/stdout - --quiet https://raw.githubusercontent.com/MesserLab/SLiM-Extras/master/installation/DebianUbuntuInstall.sh | sudo bash -s
    
mkdir -p /xdisk/mcnew/finches/dannyjackson/simulations/output/logs

SLIM_SCRIPT="/xdisk/mcnew/finches/dannyjackson/simulations/simulation_replicates.slim"



sbatch run_slim.sh -s 0.00 -r 0.95 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.50 -r 0.95 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.95 -r 0.95 -f 0.25 -m hudson

sbatch run_slim.sh -s 0.00 -r 0.85 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.50 -r 0.85 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.95 -r 0.85 -f 0.25 -m hudson

sbatch run_slim.sh -s 0.00 -r 0.75 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.50 -r 0.75 -f 0.25 -m hudson
sbatch run_slim.sh -s 0.95 -r 0.75 -f 0.25 -m hudson


module load micromamba
micromamba activate r_elgato

Rscript plot_fst_delta_tajd.R

/xdisk/mcnew/finches/dannyjackson/simulations
cat output/0.95/selection_0.00/summary_stats.tsv

cat /xdisk/mcnew/finches/dannyjackson/simulations/output/0.95/selection_0.00/replicate_1/recap_stats_run1.log | head -200

micromamba create -n recap_py \
  -c conda-forge -c bioconda \
  python=3.11 \
  msprime \
  pyslim \
  tskit \
  numpy

micromamba activate recap_py

python recap_fst_tajd.py \
  --trees output/1/selection_0.00/replicate_1/simulation_run1.trees \
  --t1_ids output/1/selection_0.00/replicate_1/sample_t1_run1.ids.txt \
  --t2_ids output/1/selection_0.00/replicate_1/sample_t2_run1.ids.txt \
  --Ne 150000 --mu 2.04e-9 --recomb 1e-8 --L 50000 \
  --model hudson \
  --seed 1 \
  --out_trees simulation_run1.recap.mut.trees \
  --verbose --suppress_time_warning


 ls output/1/selection_0.00/replicate_1/
run1.log  sample_t1_run1.ids.txt  sample_t2_run1.ids.txt  simulation_run1.trees






# sbatch run_slim.sh -s 0.00 -r 0.95
# sbatch run_slim.sh -s 0.25 -r 0.95 # running on elgato
# sbatch run_slim.ocelote.sh -s 0.5 -r 0.95 # running on puma
# sbatch run_slim.ocelote.sh -s 0.75 -r 0.95 # running on ocelote
sbatch run_slim.ocelote.sh -s 1 -r 0.95

sbatch run_slim.sh -s 0.00 -r 0.90
sbatch run_slim.sh -s 0.25 -r 0.90
sbatch run_slim.sh -s 0.5 -r 0.90
sbatch run_slim.sh -s 0.75 -r 0.90
sbatch run_slim.sh -s 1 -r 0.90

sbatch run_slim.sh -s 0.00 -r 0.85
sbatch run_slim.sh -s 0.25 -r 0.85
sbatch run_slim.sh -s 0.5 -r 0.85
sbatch run_slim.sh -s 0.75 -r 0.85
sbatch run_slim.sh -s 1 -r 0.85

# sbatch run_slim.Ne250k.sh -s 0.00 -r 0.95 # running on puma
# sbatch run_slim.Ne250k.sh -s 0.25 -r 0.95 # running on elgato
# sbatch run_slim.Ne250k.sh -s 0.5 -r 0.95 # running on ocelote

#!/usr/bin/env bash

jid=""

for r in 0.95 0.90 0.85; do
  for s in 0.25 0.5 0.75 1; do
    if [[ -z "$jid" ]]; then
      jid=$(sbatch /xdisk/mcnew/finches/dannyjackson/simulations/run_slim.sh -s "$s" -r "$r" | awk '{print $4}')
    else
      jid=$(sbatch --dependency=afterok:$jid /xdisk/mcnew/finches/dannyjackson/simulations/run_slim.sh -s "$s" -r "$r" | awk '{print $4}')
    fi
  done
done

chmod +x submit_chain.sh
./submit_chain.sh


# Create files with starting frequency of each selected allele within each selection_* dir
for d in selection_*; do
  awk '
    /picked mut/ {
      run = $2
      sub(/:$/, "", run)

      for (i=1; i<=NF; i++) {
        if ($i ~ /^freq=/) {
          freq = $i
          sub(/^freq=/, "", freq)
        }
      }
      print run "\t" freq
    }
  ' "$d"/replicate_*/run*.log 2>/dev/null | sort -u > "$d/freq.txt"

  echo "Wrote $d/freq.txt"
done


# Combine freq, FST, and Tajima's D across all selection coefficients into a single summary file

# Combine startFreq + FST + Tajima's D across all selection coefficients

echo -e "selCoef\tRun\tstartFreq\tfst\tdelta_tD" > composite_stats.tsv

for d in selection_*; do
  sel=$(echo "$d" | sed 's/selection_//')

  awk -v sel="$sel" '
    FNR==NR {
      # freq.txt: "14  0.4213"  -> key "run14"
      run = "run"$1
      freq[run] = $2
      next
    }
    FNR==NR { next }  # (not strictly necessary; keeps structure clear)

    # fst.txt: run fst
    FILENAME ~ /fst\.txt$/ {
      fst[$1] = $2
      next
    }

    # tajimaD.txt: run T1/T2 value
    {
      run = $1
      time = $2
      val = $3

      if (time == "T1") t1[run] = val
      if (time == "T2") t2[run] = val
    }

    END {
      for (run in fst) {
        if ((run in t1) && (run in t2)) {
          delta = t2[run] - t1[run]
          sf = (run in freq ? freq[run] : "NA")
          printf "%s\t%s\t%s\t%.6f\t%.6f\n", sel, run, sf, fst[run], delta
        }
      }
    }
  ' "$d/freq.txt" "$d/fst.txt" "$d/tajimaD.txt" >> composite_stats.tsv
done





module load micromamba
micromamba activate r_elgato

R

# Plot it
library(readr)
library(ggplot2)

df <- read_tsv("composite_stats.tsv", show_col_types = FALSE)

ggplot(df, aes(x = fst, y = delta_tD, color = startFreq)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 1.2, alpha = 0.85) +
  facet_wrap(~ selCoef, ncol = 3) +
  theme_classic() +
  scale_color_viridis_c(
    name = "Starting\nfrequency",
    option = "plasma",
    na.value = "grey80"
  ) +
  labs(
    x = "FST",
    y = expression(Delta*" Tajima's D"),
    title = "FST vs ΔTajima’s D across selection coefficients"
  )

ggsave(
  "composite_stats.simulated.faceted.png",
  width = 10,
  height = 5,
  dpi = 300
)

