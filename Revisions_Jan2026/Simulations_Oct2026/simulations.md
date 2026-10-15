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
```


bash ${SCRIPTDIR}/topup_selected.sh

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
# Modify plots
```
# Power curves from D1_RankPower output: one figure.
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

res <- read.delim("output/fst_rank_by_selected_windows.nnull10000.tsv")

metric      <- "prop_top_5pct"      # or prop_top_1pct / prop_top_10pct
metric_lab  <- "Power (pathway ranks in top 5%)"
null_level  <- 0.05                 # expected value for an all-neutral set
base_Ne     <- "100k"
base_decl   <- 1

# ---- 1. Do Ne and decline matter? (largest spread in power across them) ----
spread <- res %>%
  group_by(Gen, sel_s, af, X_selected) %>%
  summarise(range = diff(range(.data[[metric]])), .groups = "drop")
cat(sprintf("Max spread in %s across Ne x decline: %.2f (median %.2f)\n",
            metric, max(spread$range), median(spread$range)))

# ---- 2. Main figure at the baseline Ne / decline ----
d <- res %>%
  filter(Ne == base_Ne, decline == base_decl) %>%
  mutate(s   = factor(sel_s, levels = sort(unique(sel_s))),
         f0  = factor(paste0("f0 = ", af), levels = paste0("f0 = ", sort(unique(af)))),
         Gen = factor(paste0(Gen, " generations"),
                      levels = paste0(sort(unique(Gen)), " generations")))

s_cols <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")
s_cols <- s_cols[seq_len(nlevels(d$s))]

p <- ggplot(d, aes(X_selected, .data[[metric]], colour = s, group = s)) +
  geom_hline(yintercept = null_level, linetype = "dashed", colour = "#8a8984", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) +
  facet_grid(Gen ~ f0) +
  scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8, 10)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Selected windows in a 10-window pathway",
       y = metric_lab,
       title = "Power to detect a selected pathway",
       subtitle = sprintf("Ne = %s, no decline; dashed line = expectation for an all-neutral pathway",
                          base_Ne)) +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        strip.text = element_text(face = "bold", colour = "#3d3d3a"),
        axis.text = element_text(colour = "#8a8984"),
        legend.position = "right")

dir.create("output/power_curves", showWarnings = FALSE, recursive = TRUE)
ggsave("output/power_curves/power_curves_main.png", p, width = 13, height = 6, dpi = 300, bg = "white")

# ---- 3. All conditions, grouped Ne > generations > decline ----
ne_num <- function(x) ifelse(grepl("k$", x), as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

d_all <- res %>%
  mutate(
    Ne_num  = ne_num(Ne),
    s       = factor(sel_s, levels = sort(unique(sel_s))),
    f0      = factor(paste0("f0 = ", af), levels = paste0("f0 = ", sort(unique(af)))),
    Ne_lab  = factor(paste0("Ne ", Ne), levels = paste0("Ne ", unique(Ne[order(Ne_num)]))),
    Gen_lab = factor(paste0(Gen, " gen"), levels = paste0(sort(unique(Gen)), " gen")),
    # no decline first, then declining
    Dec_lab = factor(ifelse(decline == 1, "no decline", paste0("decline ", decline)),
                     levels = c("no decline",
                                paste0("decline ", sort(unique(decline[decline != 1]), decreasing = TRUE))))
  )

curve_layers <- function(base_size, lw) list(
  geom_hline(yintercept = null_level, linetype = "dashed", colour = "#8a8984", linewidth = 0.3),
  geom_line(linewidth = lw),
  scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)"),
  scale_x_continuous(breaks = c(0, 5, 10)),
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)),
  labs(x = "Selected windows in a 10-window pathway", y = metric_lab),
  theme_minimal(base_size = base_size),
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        strip.text = element_text(face = "bold", colour = "#3d3d3a"),
        strip.text.y = element_text(angle = 0, hjust = 0))
)

# 3a. One figure: rows nested Ne > Gen > decline, columns f0
n_rows <- nlevels(d_all$Ne_lab) * nlevels(d_all$Gen_lab) * nlevels(d_all$Dec_lab)
p_all <- ggplot(d_all, aes(X_selected, .data[[metric]], colour = s, group = s)) +
  curve_layers(9, 0.5) +
  facet_grid(Ne_lab + Gen_lab + Dec_lab ~ f0)
ggsave("output/power_curves/power_curves_all_conditions.png", p_all,
       width = 13, height = 1.4 * n_rows + 1.5, dpi = 300, bg = "white", limitsize = FALSE)

# 3b. One file per combination, filed Ne_*/Gen_*/, decline in the file name
combos <- d_all %>% distinct(Ne, Ne_num, Gen, decline, Ne_lab, Gen_lab, Dec_lab) %>%
  arrange(Ne_num, Gen, desc(decline))

for (i in seq_len(nrow(combos))) {
  cb <- combos[i, ]
  di <- d_all %>% filter(Ne == cb$Ne, Gen == cb$Gen, decline == cb$decline)
  pi <- ggplot(di, aes(X_selected, .data[[metric]], colour = s, group = s)) +
    curve_layers(11, 0.7) +
    geom_point(size = 1.1) +
    facet_grid(. ~ f0) +
    labs(title = sprintf("Ne %s | %s generations | %s", cb$Ne, cb$Gen, cb$Dec_lab))
  outdir <- file.path("output/power_curves", paste0("Ne_", cb$Ne), paste0("Gen_", cb$Gen))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(outdir, sprintf("power_curves.decline_%.2f.png", cb$decline)), pi,
         width = 13, height = 3.6, dpi = 300, bg = "white")
}

cat("Saved output/power_curves/power_curves_main.png,",
    "power_curves_all_conditions.png and per-condition plots in output/power_curves/Ne_*/Gen_*/\n")
```