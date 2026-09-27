# Modeling sensitivity of our sampling scheme and analytical approach

cd /xdisk/mcnew/finches/dannyjackson/simulations/power_analysis_ttp

## Model neutral evolution
```
# 10000 50000 100000

F0=0.05
S=0.0
L=10000
MODEL="hudson"

NE_VALUES=(10000)
DECLINES=(1.00)

# Generations between t1 and t2
GEN_TIMES_T1_T2=(15)

# Generations between t2 and t3
GEN_TIMES_T2_T3=(5)

for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G12 in "${GEN_TIMES_T1_T2[@]}"; do
      for G23 in "${GEN_TIMES_T2_T3[@]}"; do

        # The SLiM script expects T3_OFFSET measured from T1,
        # so convert the T2->T3 interval to cumulative T1->T3 time.
        G13=$((G12 + G23))

        echo "Submitting: Ne=${NE} decline=${D} T1-T2=${G12} T2-T3=${G23} T1-T3=${G13}"

        sbatch run_slim.power_analysis.ttp.sh \
          -s "${S}" \
          -r "${D}" \
          -f "${F0}" \
          -m "${MODEL}" \
          -N "${NE}" \
          -o "${G12}" \
          -O "${G13}" \
          -L "${L}" \
          -n 9900

      done
    done
  done
done

```
## Model positive selection
```
MODEL="hudson"
F0S=(0.05 0.5)
NES=(10000 50000 100000)
SELS=(0.01 0.05 0.10 0.20 0.30 0.40 0.50)
DECLINES=(1.00 0.80)

# Generations between t1 and t2
GEN_TIMES_T1_T2=(5 15 30)

# Generations between t2 and t3
GEN_TIMES_T2_T3=(5)

LENS=(10000)

for F0 in "${F0S[@]}"; do
  for NE in "${NES[@]}"; do
    for S in "${SELS[@]}"; do
      for D in "${DECLINES[@]}"; do
        for G12 in "${GEN_TIMES_T1_T2[@]}"; do
          for G23 in "${GEN_TIMES_T2_T3[@]}"; do
            for L in "${LENS[@]}"; do

              # T3_OFFSET is measured from T1
              G13=$((G12 + G23))

              echo "Submitting: f0=${F0} Ne=${NE} s=${S} decline=${D} T1-T2=${G12} T2-T3=${G23} T1-T3=${G13} L=${L}"

              sbatch run_slim.power_analysis.sh \
                -s "${S}" \
                -r "${D}" \
                -f "${F0}" \
                -m "${MODEL}" \
                -N "${NE}" \
                -o "${G12}" \
                -O "${G13}" \
                -L "${L}" \
                -n 100

              sleep 0.1

            done
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

# Rscript plot.power_analysis.r
# Rscript plot.power_analysis.f0_join.r
Rscript plot.power_analysis.linkf0.r
```
library(dplyr)
library(ggplot2)
library(tidyr)

df <- read.csv("plots_fst_deltaTajD_neutral_vs_selected/power_comparison_wide.tsv", sep = '\t')

# Convert to long format for statistics
df_long <- df %>%
  pivot_longer(
    cols = c(power_composite, power_fst, power_taj),
    names_to = "statistic",
    values_to = "power"
  )

df_long <- df_long %>%
  mutate(
    Ne = factor(Ne),
    decline = factor(decline),
    Gen_time = factor(Gen_time),
    sel_s = factor(sel_s),
    f0_selected = factor(f0_selected)
  )

# Unique combinations
combos <- df_long %>%
  distinct(f0_selected, statistic)

for(i in seq_len(nrow(combos))){

  f0_val <- combos$f0_selected[i]
  stat_val <- combos$statistic[i]

  df_plot <- df_long %>%
    filter(
      f0_selected == f0_val,
      statistic == stat_val
    )

  p <- ggplot(df_plot, aes(x = sel_s, y = Gen_time, fill = power)) +
    geom_tile(color = "white") +
    facet_grid(decline ~ Ne, labeller = label_both) +
    scale_fill_viridis_c(
      option = "magma",
      name = "Power",
      limits = c(0,1),
      breaks = seq(0,1,by=0.2),
      oob = scales::squish
    ) +
    labs(
      x = "Selection coefficient (s)",
      y = "Generations between samples",
      title = paste("Power:", stat_val, "| f0 =", f0_val)
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90")
    )

  ggsave(
    paste0(
      "plots_fst_deltaTajD_neutral_vs_selected/power.",
      stat_val,
      ".f0_", f0_val,
      ".png"
    ),
    plot = p,
    width = 10,
    height = 3
  )

}
```