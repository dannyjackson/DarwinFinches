#!/usr/bin/env Rscript
# ============================================================
# Publication figure: power of four approaches (2 x 2 grid)
#
#                 gene level (outliers)          pathway level (GSEA)
#   one species   A  single-species outliers     B  single-species GSEA
#   three species C  three-species overlap       D  cross-species GSEA
#                                                   (blood coagulation)
#
# Every panel: x = selection coefficient (log), y = power, one line per
# starting frequency f0, dashed line = expectation with no selection.
# One figure per generation time. Ne and decline fixed (see settings).
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/plot_power_summary_figure.R
#
# Inputs (paths below):
#   A  output/fst_outlier_power.nnull9900.nsel100.tsv
#   B  output/fst_rank_by_selected_windows.nnull10000.tsv
#   C  output/parallel_overlap.permutations.tsv
#   D  gsea_empirical_power/summary.tsv
#
# Writes: output/figure_power_summary/power_summary.Gen_<g>.png / .pdf
#         output/figure_power_summary/power_summary.data.tsv  (values plotted)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---------------- settings ----------------
file_outlier <- "output/fst_outlier_power.nnull9900.nsel100.tsv"
file_gsea1   <- "output/fst_rank_by_selected_windows.nnull10000.tsv"
file_overlap <- "output/parallel_overlap.permutations.tsv"
file_gsea3   <- "gsea_empirical_power/summary.tsv"

out_dir <- "output/figure_power_summary"

base_Ne      <- "100k"
base_decline <- 1

# GSEA threshold: "5" (top 5%) or "1" (top 1%) for panels B and D
gsea_top <- "5"

# Panel B: number of selected windows in the 10-window pathway
gsea1_X <- 10

# Panel C: "recall" = share of the 30 parallel genes found in all three top lists
#          "any"    = P(at least one parallel gene found in all three top lists)
overlap_metric <- "recall"
overlap_n_parallel <- 30

text_size <- 12

# light (low f0) -> dark (high f0)
f0_ramp <- c("#9ec5f4", "#5598e7", "#2a78d6", "#184f95", "#0d366b")

# ---------------- helpers ----------------
norm_keys <- function(d, f0_col) {
  d %>%
    mutate(Ne = as.character(Ne),
           Gen = as.numeric(Gen),
           decline = as.numeric(decline),
           sel_s = as.numeric(sel_s),
           f0 = as.numeric(.data[[f0_col]]))
}

base_filter <- function(d) {
  d %>% filter(Ne == base_Ne, abs(decline - base_decline) < 1e-9, sel_s > 0)
}

panel_lab <- c(
  A = "A  One species: outlier windows",
  B = "B  One species: pathway enrichment",
  C = "C  Three species: outlier overlap",
  D = "D  Three species: pathway enrichment (blood coagulation)"
)

# ---------------- A: single-species outliers ----------------
A <- read.delim(file_outlier) %>% norm_keys("af") %>% base_filter() %>%
  transmute(panel = "A", Gen, sel_s, f0, power = mean_power)
null_A <- 0.01

# ---------------- B: single-species GSEA ----------------
col_B <- paste0("prop_top_", gsea_top, "pct")
B <- read.delim(file_gsea1) %>% norm_keys("af") %>% base_filter() %>%
  filter(X_selected == gsea1_X) %>%
  transmute(panel = "B", Gen, sel_s, f0, power = .data[[col_B]])
null_B <- as.numeric(gsea_top) / 100

# ---------------- C: three-species overlap ----------------
C <- read.delim(file_overlap) %>% norm_keys("af") %>% base_filter() %>%
  filter(n_parallel == overlap_n_parallel) %>%
  group_by(Gen, sel_s, f0) %>%
  summarise(power = if (overlap_metric == "recall") mean(parallel_in_all_tops)
                    else mean(overlap3_parallel > 0),
            .groups = "drop") %>%
  mutate(panel = "C")
null_C <- 0   # chance of a gene in all three top-1% lists: 0.01^3

# ---------------- D: cross-species GSEA, blood coagulation ----------------
col_D <- paste0("power_cross_top", gsea_top)
D_all <- read.delim(file_gsea3) %>% norm_keys("f0") %>%
  filter(Ne == base_Ne, abs(decline - base_decline) < 1e-9)
D <- D_all %>% filter(sel_s > 0) %>%
  transmute(panel = "D", Gen, sel_s, f0, power = .data[[col_D]])
null_D <- as.numeric(gsea_top) / 100   # observed s = 0 rows printed below for comparison

if (any(D_all$sel_s == 0)) {
  obs <- D_all %>% filter(sel_s == 0) %>% select(Gen, all_of(col_D))
  cat("Panel D observed null (s = 0):\n"); print(obs)
}

# ---------------- combine ----------------
dat <- bind_rows(A, B, C, D) %>%
  mutate(panel_label = factor(panel_lab[panel], levels = panel_lab),
         f0_lab = factor(f0, levels = sort(unique(f0))))

nulls <- tibble(panel = c("A", "B", "C", "D"),
                null = c(null_A, null_B, null_C, null_D)) %>%
  mutate(panel_label = factor(panel_lab[panel], levels = panel_lab))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(panel, Gen, sel_s, f0, power) %>% arrange(panel, Gen, f0, sel_s),
            file.path(out_dir, "power_summary.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Report which s values each panel has
cat("\ns values per panel:\n")
print(dat %>% group_by(panel) %>%
        summarise(s = paste(sort(unique(sel_s)), collapse = ", "), .groups = "drop"))

f0_cols <- setNames(rep_len(f0_ramp, nlevels(dat$f0_lab)), levels(dat$f0_lab))
s_breaks <- sort(unique(dat$sel_s))

# Panel subtitles: what "power" means in each
defs <- tibble(
  panel = c("A", "B", "C", "D"),
  def = c(
    "Share of selected windows in the top 1% (100 of 10,000)",
    sprintf("P(pathway with %d/10 selected windows ranks in top %s%%)", gsea1_X, gsea_top),
    if (overlap_metric == "recall")
      sprintf("Share of %d parallel genes in all three top-1%% lists", overlap_n_parallel)
    else sprintf("P(>=1 of %d parallel genes in all three top-1%% lists)", overlap_n_parallel),
    sprintf("P(blood coagulation ranks in top %s%% by mean NES)", gsea_top)
  )
) %>% mutate(panel_label = factor(panel_lab[panel], levels = panel_lab))

# ---------------- plot ----------------
for (g in sort(unique(dat$Gen))) {

  dg <- dat %>% filter(Gen == g)
  if (nrow(dg) == 0) next

  p <- ggplot(dg, aes(sel_s, power, colour = f0_lab, group = f0_lab)) +
    geom_hline(data = nulls, aes(yintercept = null),
               linetype = "dashed", colour = "#8a8984", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    geom_text(data = defs, aes(x = min(s_breaks), y = 1, label = def),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = text_size * 0.26,
              colour = "#6b6a66") +
    facet_wrap(~ panel_label, ncol = 2) +
    scale_colour_manual(values = f0_cols, name = "Starting allele\nfrequency (f0)") +
    scale_x_log10(breaks = s_breaks, labels = as.character(s_breaks)) +
    scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
    labs(x = "Selection coefficient (s, log scale)", y = "Power",
         title = sprintf("Power to detect selection after %s generations", g),
         subtitle = sprintf("Ne = %s, %s. Dashed line: expectation with no selection.",
                            base_Ne, ifelse(base_decline == 1, "constant size",
                                            paste0("decline ", base_decline, " per generation")))) +
    theme_minimal(base_size = text_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
      panel.spacing    = unit(1.2, "lines"),
      strip.text       = element_text(face = "bold", colour = "#3d3d3a", hjust = 0,
                                      size = text_size),
      axis.text        = element_text(colour = "#6b6a66"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      plot.title       = element_text(face = "bold"),
      legend.position  = "right"
    )

  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, sprintf("power_summary.Gen_%s.%s", g, ext)), p,
           width = 11, height = 8, dpi = 300, bg = "white")
  }
  message("Saved: ", file.path(out_dir, sprintf("power_summary.Gen_%s.png", g)))
}