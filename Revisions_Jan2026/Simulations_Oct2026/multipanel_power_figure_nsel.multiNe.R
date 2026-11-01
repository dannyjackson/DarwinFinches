#!/usr/bin/env Rscript
# ============================================================
# Alternate publication figure: power vs. the NUMBER OF SELECTED GENES (N)
# for all four approaches (2 x 2, same arrangement as multipanel_power_figure.R)
#
#                     one species                    three species (parallel)
#   outliers (top)    A  selected windows in genome  B  parallel genes in genome
#   GSEA (bottom)     C  selected windows in a       D  selected genes in blood
#                        10-window pathway              coagulation
#
# Every panel: x = N, y = power, one line per Ne; one f0 and one s per figure
# (settings below).
# The gene-level tests (A, B) score each gene on its own, so per-gene power is
# expected to stay flat as N grows; the pathway tests (C, D) pool genes, so
# power should rise with N. That contrast is the point of this figure.
#
# Inputs (run these first):
#   A  output/fst_outlier_power_nsel.tsv        outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv parallel_outlier_overlap.R with
#                                               PARALLEL_VALUES <- c(1, 2, 5, 10, 20, 30, 50, 100)
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv          plot_gsea_parallel_nsel.R
#
# Run from the project directory:
#   Rscript ${SCRIPTDIR}/multipanel_power_figure_nsel.R
# Writes: output/figure_power_summary_nsel/power_vs_n.byNe.Gen_<g>.f0_<f0>.s_<s>.png / .pdf
#         output/figure_power_summary_nsel/power_vs_n.byNe.data.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)   # install.packages("patchwork") if missing
})

# ---------------- settings ----------------
file_outlier <- "output/fst_outlier_power_nsel.tsv"
file_overlap <- "output/parallel_overlap.permutations.tsv"
file_gsea1   <- "output/fst_rank_by_selected_windows.nnull10000.tsv"
file_gsea3   <- "gsea_empirical_power_nsel/summary_nsel.tsv"
out_dir      <- "output/figure_power_summary_nsel"

Ne_values    <- c("25k", "50k", "75k")   # one line per Ne in every panel
base_decline <- 1
f0_fixed     <- 0.2          # one starting frequency for every panel
s_values     <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5)   # one figure per s (per Gen); NULL = every s
gsea_top     <- "5"          # GSEA cutoff for C and D: "5" or "1" (top %)

text_size <- 12

# one distinct hue per Ne (small -> large), colour-blind safe (Okabe-Ito)
Ne_pal <- c("25k" = "#E69F00", "50k" = "#0072B2", "75k" = "#009E73")

# ---------------- helpers ----------------
prep <- function(d, f0_col) {
  d %>% mutate(Ne = as.character(Ne), Gen = as.numeric(Gen),
               decline = as.numeric(decline), sel_s = as.numeric(sel_s),
               f0 = as.numeric(.data[[f0_col]])) %>%
    filter(Ne %in% Ne_values, abs(decline - base_decline) < 1e-9,
           abs(f0 - f0_fixed) < 1e-9, sel_s > 0)
}

# ---------------- A: single-species outliers vs. N selected windows ----------------
A <- read.delim(file_outlier) %>% prep("af") %>%
  transmute(panel = "A", Ne, Gen, sel_s, N = n_selected, power = recall)

# ---------------- B: three-species overlap vs. N parallel genes ----------------
B <- read.delim(file_overlap) %>% prep("af") %>%
  filter(n_parallel > 0) %>%
  group_by(Ne, Gen, sel_s, N = n_parallel) %>%
  summarise(power = mean(parallel_in_all_tops), .groups = "drop") %>%
  mutate(panel = "B")
if (n_distinct(B$N) < 2)
  warning("Panel B has only N = ", paste(unique(B$N), collapse = ", "),
          ". Rerun parallel_outlier_overlap.R with several PARALLEL_VALUES.")

# ---------------- C: single-species GSEA vs. X selected windows ----------------
C <- read.delim(file_gsea1) %>% prep("af") %>%
  transmute(panel = "C", Ne, Gen, sel_s, N = X_selected,
            power = .data[[paste0("prop_top_", gsea_top, "pct")]])

# ---------------- D: cross-species GSEA vs. N selected genes ----------------
D_raw <- read.delim(file_gsea3)
D <- D_raw %>% prep("f0") %>%
  transmute(panel = "D", Ne, Gen, sel_s, N = n_selected,
            power = .data[[paste0("power_cross_top", gsea_top)]])
n_term <- max(D_raw$n_selected)

dat <- bind_rows(A, B, C, D) %>% mutate(Ne_lab = factor(Ne, levels = Ne_values))
if (!is.null(s_values)) dat <- dat %>% filter(round(sel_s, 6) %in% round(s_values, 6))

# Flag any panel x Ne with no rows (e.g. an upstream script not yet run for that Ne)
miss <- expand.grid(panel = c("A", "B", "C", "D"), Ne = Ne_values, stringsAsFactors = FALSE) %>%
  anti_join(distinct(dat, panel, Ne), by = c("panel", "Ne"))
if (nrow(miss) > 0) warning("No data for: ",
                            paste(miss$panel, miss$Ne, sep = "/", collapse = ", "))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(-Ne_lab) %>% arrange(panel, Gen, sel_s, Ne, N),
            file.path(out_dir, "power_vs_n.byNe.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("N values per panel:\n")
print(dat %>% group_by(panel) %>%
        summarise(N = paste(sort(unique(N)), collapse = ", "),
                  s = paste(sort(unique(sel_s)), collapse = ", "), .groups = "drop"))

Ne_cols <- Ne_pal[Ne_values]

# ---------------- panel specs ----------------
specs <- list(
  A = list(title = "A  Outliers: one species",
           def   = "Share of selected windows in the top 1% of 10,000",
           xlab  = "Selected windows in the genome (N, log scale)",
           null  = 0.01, log = TRUE),
  B = list(title = "B  Outliers: parallel across three species",
           def   = "Share of parallel genes in all three top-1% lists",
           xlab  = "Genes selected in all three species (N, log scale)",
           null  = 0, log = TRUE),
  C = list(title = "C  Pathway enrichment: one species",
           def   = sprintf("P(10-window pathway ranks in top %s%%)", gsea_top),
           xlab  = "Selected windows in the pathway (N, of 10)",
           null  = as.numeric(gsea_top) / 100, log = FALSE),
  D = list(title = "D  Pathway enrichment: parallel across three species",
           def   = sprintf("P(blood coagulation in top %s%% by mean NES)", gsea_top),
           xlab  = sprintf("Selected blood-coagulation genes (N, of %d)", n_term),
           null  = as.numeric(gsea_top) / 100, log = FALSE)
)

make_panel <- function(d, sp) {
  xs <- sort(unique(d$N))
  p <- ggplot(d, aes(N, power, colour = Ne_lab, group = Ne_lab)) +
    geom_hline(yintercept = sp$null, linetype = "dashed", colour = "#8a8984", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = Ne_cols, name = "Effective\npopulation\nsize (Ne)", drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = sp$title, subtitle = sp$def, x = sp$xlab, y = "Power") +
    theme_minimal(base_size = text_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
          axis.text = element_text(colour = "#6b6a66"),
          plot.title = element_text(face = "bold", size = text_size),
          plot.subtitle = element_text(colour = "#6b6a66", size = text_size * 0.85))
  if (sp$log) p + scale_x_log10(breaks = xs, labels = as.character(xs))
  else        p + scale_x_continuous(breaks = xs)
}

for (g in sort(unique(dat$Gen))) {
 for (s_now in sort(unique(dat$sel_s))) {
  dg <- dat %>% filter(Gen == g, abs(sel_s - s_now) < 1e-9)
  if (nrow(dg) == 0) next
  panels <- lapply(c("A", "B", "C", "D"), function(k) {
    dk <- dg %>% filter(panel == k)
    if (nrow(dk) == 0) {
      ggplot() + labs(title = specs[[k]]$title, subtitle = "No data for this setting") +
        theme_void(base_size = text_size)
    } else make_panel(dk, specs[[k]])
  })

  fig <- wrap_plots(panels, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Power vs. number of selected genes after %s generations", g),
      subtitle = sprintf(paste0("s = %g, f0 = %g, %s. Dashed lines: expectation with no selection. ",
                                "Gene-level tests (A, B) score genes one at a time; ",
                                "pathway tests (C, D) pool them."),
                         s_now, f0_fixed, ifelse(base_decline == 1, "constant size",
                                                 paste0("decline ", base_decline))),
      theme = theme(plot.title = element_text(face = "bold", size = text_size + 2),
                    plot.subtitle = element_text(colour = "#6b6a66", size = text_size * 0.9)))

  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, sprintf("power_vs_n.byNe.Gen_%s.f0_%g.s_%g.%s", g, f0_fixed, s_now, ext)),
           fig, width = 12, height = 9, dpi = 300, bg = "white")
  }
  message("Saved: ", file.path(out_dir, sprintf("power_vs_n.byNe.Gen_%s.f0_%g.s_%g.png",
                                                g, f0_fixed, s_now)))
 }
}
