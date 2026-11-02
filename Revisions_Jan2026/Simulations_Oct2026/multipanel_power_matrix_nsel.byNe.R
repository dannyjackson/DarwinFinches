#!/usr/bin/env Rscript
# ============================================================
# Power matrix per Ne: all four approaches x every starting frequency
#
# One figure per Ne x decline x generation time.
#
#   rows    = approach
#               A  outliers: one species          (x = selected windows in genome)
#               B  outliers: three-species overlap (x = parallel genes in genome)
#               C  GSEA: one species              (x = selected windows in 10-window pathway)
#               D  GSEA: three species            (x = selected blood-coagulation genes)
#   columns = starting frequency f0
#   x       = number selected (N); log scale for A and B
#   y       = power
#   lines   = selection coefficient s (light = weak, dark = strong)
#   dashed  = expectation with no selection
#
# Same inputs and power definitions as multipanel_power_figure_nsel.multiNe.R:
#   A  output/fst_outlier_power_nsel.tsv                    outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv             parallel_outlier_overlap.R with
#                                                           PARALLEL_VALUES <- c(1, 2, 5, 10, 20, 30, 50, 100)
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv           plot_gsea_parallel_nsel.R
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/multipanel_power_matrix_nsel.byNe.R
#
# Writes: output/figure_power_matrix_nsel/
#           power_matrix.Ne_<Ne>.decline_<d>.Gen_<g>.png / .pdf
#           power_matrix.data.tsv   (values plotted)
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
out_dir      <- "output/figure_power_matrix_nsel"

# NULL = every value present in the inputs; or e.g. c("25k", "50k", "75k")
Ne_values      <- NULL
decline_values <- NULL   # e.g. 1 for constant size only
f0_values      <- NULL   # columns; e.g. c(0.05, 0.2, 0.6)
s_values       <- NULL   # lines;   e.g. c(0.05, 0.1, 0.2)

gsea_top <- "5"          # GSEA cutoff for C and D: "5" or "1" (top %)

text_size <- 11

# light = weak selection, dark = strong selection (up to 9 s values)
s_ramp <- c("#b7d3f6", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#1c5cab", "#104281", "#0d366b")

# ---------------- helpers ----------------
ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000,
                             as.numeric(x))

prep <- function(d, f0_col) {
  d %>%
    mutate(Ne = as.character(Ne), Gen = as.numeric(Gen),
           decline = as.numeric(decline), sel_s = as.numeric(sel_s),
           f0 = as.numeric(.data[[f0_col]])) %>%
    filter(sel_s > 0)
}

# ---------------- A: single-species outliers vs. N selected windows ----------------
A <- read.delim(file_outlier) %>% prep("af") %>%
  transmute(panel = "A", Ne, decline, Gen, sel_s, f0, N = n_selected, power = recall)

# ---------------- B: three-species overlap vs. N parallel genes ----------------
B <- read.delim(file_overlap) %>% prep("af") %>%
  filter(n_parallel > 0, n_parallel != 30) %>%
  group_by(Ne, decline, Gen, sel_s, f0, N = n_parallel) %>%
  summarise(power = mean(parallel_in_all_tops), .groups = "drop") %>%
  mutate(panel = "B")

# ---------------- C: single-species GSEA vs. X selected windows ----------------
C <- read.delim(file_gsea1) %>% prep("af") %>%
  transmute(panel = "C", Ne, decline, Gen, sel_s, f0, N = X_selected,
            power = .data[[paste0("prop_top_", gsea_top, "pct")]])

# ---------------- D: cross-species GSEA vs. N selected genes ----------------
D_raw <- read.delim(file_gsea3)
D <- D_raw %>% prep("f0") %>%
  transmute(panel = "D", Ne, decline, Gen, sel_s, f0, N = n_selected,
            power = .data[[paste0("power_cross_top", gsea_top)]])
n_term <- max(D_raw$n_selected)

# ---------------- combine and filter ----------------
dat <- bind_rows(A, B, C, D)

if (!is.null(Ne_values))      dat <- dat %>% filter(Ne %in% Ne_values)
if (!is.null(decline_values)) dat <- dat %>% filter(round(decline, 6) %in% round(decline_values, 6))
if (!is.null(f0_values))      dat <- dat %>% filter(round(f0, 6) %in% round(f0_values, 6))
if (!is.null(s_values))       dat <- dat %>% filter(round(sel_s, 6) %in% round(s_values, 6))

if (nrow(dat) == 0) stop("No rows left after filtering; check the settings.")

# Fixed levels so every row and every figure use the same columns and colours
f0_levels <- sort(unique(dat$f0))
s_levels  <- sort(unique(dat$sel_s))
if (length(s_levels) > length(s_ramp))
  stop("More s values (", length(s_levels), ") than colours in s_ramp (", length(s_ramp), ")")
s_cols <- setNames(s_ramp[round(seq(1, length(s_ramp), length.out = length(s_levels)))],
                   as.character(s_levels))

dat <- dat %>%
  mutate(f0_lab = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", f0_levels)),
         s_lab  = factor(as.character(sel_s), levels = as.character(s_levels)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(panel, Ne, decline, Gen, f0, sel_s, N, power) %>%
              arrange(Ne, decline, Gen, panel, f0, sel_s, N),
            file.path(out_dir, "power_matrix.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage report: which Ne / decline / Gen each approach has
cat("Coverage per approach:\n")
print(dat %>% group_by(panel) %>%
        summarise(Ne      = paste(unique(Ne[order(ne_num(Ne))]), collapse = ", "),
                  decline = paste(sort(unique(decline)), collapse = ", "),
                  Gen     = paste(sort(unique(Gen)), collapse = ", "),
                  f0      = paste(sort(unique(f0)), collapse = ", "),
                  s       = paste(sort(unique(sel_s)), collapse = ", "),
                  .groups = "drop"),
      width = Inf)

# ---------------- row specs ----------------
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

theme_matrix <- function() {
  theme_minimal(base_size = text_size) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
          panel.spacing    = unit(0.9, "lines"),
          strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
          axis.text        = element_text(colour = "#6b6a66"),
          plot.title       = element_text(face = "bold", size = text_size),
          plot.subtitle    = element_text(colour = "#6b6a66", size = text_size * 0.85))
}

# One row of the matrix: one approach, columns = f0
make_row <- function(d, sp) {
  xs <- sort(unique(d$N))
  if (sp$log) xs <- xs[xs > 0]
  p <- ggplot(d, aes(N, power, colour = s_lab, group = s_lab)) +
    geom_hline(yintercept = sp$null, linetype = "dashed", colour = "#8a8984", linewidth = 0.4) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.2) +
    facet_wrap(~ f0_lab, nrow = 1, drop = FALSE) +
    scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(title = sp$title, subtitle = sp$def, x = sp$xlab, y = "Power") +
    theme_matrix()
  if (sp$log) p + scale_x_log10(breaks = xs, labels = as.character(xs))
  else        p + scale_x_continuous(breaks = xs)
}

empty_row <- function(sp, why) {
  ggplot() + labs(title = sp$title, subtitle = why) + theme_void(base_size = text_size) +
    theme(plot.title = element_text(face = "bold", size = text_size),
          plot.subtitle = element_text(colour = "#6b6a66", size = text_size * 0.85))
}

# ---------------- plot: one file per Ne x decline x Gen ----------------
combos <- dat %>% distinct(Ne, decline, Gen) %>% arrange(ne_num(Ne), desc(decline), Gen)

for (i in seq_len(nrow(combos))) {
  ne <- combos$Ne[i]; dec <- combos$decline[i]; g <- combos$Gen[i]

  dg <- dat %>% filter(Ne == ne, abs(decline - dec) < 1e-9, Gen == g)

  rows <- lapply(c("A", "B", "C", "D"), function(k) {
    dk <- dg %>% filter(panel == k)
    if (nrow(dk) == 0) empty_row(specs[[k]], "No data for this Ne / decline / generation time")
    else make_row(dk, specs[[k]])
  })

  fig <- wrap_plots(rows, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Power vs. number selected: Ne = %s, %s, %s generations",
                      ne, ifelse(dec == 1, "constant size", paste0("decline ", dec)), g),
      subtitle = paste0("Columns: starting frequency f0. Lines: selection coefficient s. ",
                        "Dashed lines: expectation with no selection."),
      theme = theme(plot.title = element_text(face = "bold", size = text_size + 3),
                    plot.subtitle = element_text(colour = "#6b6a66", size = text_size)))

  w <- max(8, 2.6 * length(f0_levels) + 2)
  stem <- sprintf("power_matrix.Ne_%s.decline_%.2f.Gen_%s", ne, dec, g)
  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, paste0(stem, ".", ext)), fig,
           width = w, height = 13, dpi = 300, bg = "white")
  }
  message("Saved: ", file.path(out_dir, paste0(stem, ".png")))
}

message("\nDone. Figures written to: ", out_dir)
