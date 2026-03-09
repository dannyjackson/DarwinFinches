# Modeling sensitivity of our sampling scheme and analytical approach

cd /xdisk/mcnew/finches/dannyjackson/simulations/power_analysis/

## Model neutral evolution
```
# 10000 50000 100000
F0=0.5
S=0.0
L=10000
MODEL="hudson"


NE_VALUES=(10000)
DECLINES=(1.00)
GEN_TIMES=(5 15 30)

for NE in "${NE_VALUES[@]}"; do
  for D in "${DECLINES[@]}"; do
    for G in "${GEN_TIMES[@]}"; do

      echo "Submitting: Ne=${NE} decline=${D} gen_time=${G}"

      sbatch run_slim.power_analysis.sh \
        -s "${S}" \
        -r "${D}" \
        -f "${F0}" \
        -m "${MODEL}" \
        -N "${NE}" \
        -o "${G}" \
        -L "${L}" \
        -n 9900

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
GEN_TIMES=(5 15 30)
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