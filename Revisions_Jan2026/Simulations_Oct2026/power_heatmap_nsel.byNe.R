#!/usr/bin/env Rscript
# ============================================================
# Power heatmaps per Ne: all four approaches x generation time x f0
#
# One figure per Ne x decline.
#
#   row blocks = approach (each block has its own N axis)
#                  A  outliers: one species           (N = selected windows in genome)
#                  B  outliers: three-species overlap (N = parallel genes in genome)
#                  C  GSEA: one species               (N = selected windows in 10-window pathway)
#                  D  GSEA: three species             (N = selected blood-coagulation genes)
#   rows       = generation time, within each block
#   columns    = starting frequency f0
#   tile       = x: number selected (N), y: selection coefficient s, fill: power
#   outline    = tiles with power >= power_goal (default 0.8)
#
# The empirical N is unknown, so every N in each sweep is kept.
# Colour scale is fixed at 0-1 so colours mean the same thing in every file.
# Each block's subtitle gives the power expected with no selection.
#
# Same inputs and power definitions as multipanel_power_figure_nsel.multiNe.R:
#   A  output/fst_outlier_power_nsel.tsv                    outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv             parallel_outlier_overlap.R with
#                                                           PARALLEL_VALUES <- c(1, 2, 5, 10, 20, 30, 50, 100)
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv           plot_gsea_parallel_nsel.R
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/power_heatmap_nsel.byNe.R
#
# Writes: output/figure_power_heatmap_nsel/
#           power_heatmap.Ne_<Ne>.decline_<d>.png / .pdf
#           power_heatmap.data.tsv   (values plotted)
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
out_dir      <- "output/figure_power_heatmap_nsel"

# NULL = every value present in the inputs; or e.g. c("25k", "50k", "75k")
Ne_values      <- NULL
decline_values <- NULL   # e.g. 1 for constant size only
f0_values      <- NULL   # columns; e.g. c(0.05, 0.2, 0.6)
s_values       <- NULL   # tile rows; e.g. c(0.05, 0.1, 0.2)

gsea_top    <- "5"       # GSEA cutoff for C and D: "5" or "1" (top %)
power_goal  <- 0.8       # outline tiles at or above this power; NA = no outline
show_values <- TRUE      # print the power value in each tile

text_size <- 11

# light = low power, dark = high power
fill_ramp <- c("#f7f9fc", "#b7d3f6", "#6da7ec", "#2a78d6", "#104281", "#0d366b")

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
  filter(n_parallel > 0) %>%
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

# Fixed levels so every block and every figure use the same rows, columns and order
f0_levels  <- sort(unique(dat$f0))
s_levels   <- sort(unique(dat$sel_s))
gen_levels <- sort(unique(dat$Gen))

dat <- dat %>%
  mutate(f0_lab  = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", f0_levels)),
         s_lab   = factor(as.character(sel_s), levels = as.character(s_levels)),
         Gen_lab = factor(paste0(Gen, " gen"), levels = paste0(gen_levels, " gen")))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(panel, Ne, decline, Gen, f0, sel_s, N, power) %>%
              arrange(Ne, decline, panel, Gen, f0, sel_s, N),
            file.path(out_dir, "power_heatmap.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage report: which values each approach has
cat("Coverage per approach:\n")
print(dat %>% group_by(panel) %>%
        summarise(Ne      = paste(unique(Ne[order(ne_num(Ne))]), collapse = ", "),
                  decline = paste(sort(unique(decline)), collapse = ", "),
                  Gen     = paste(sort(unique(Gen)), collapse = ", "),
                  f0      = paste(sort(unique(f0)), collapse = ", "),
                  s       = paste(sort(unique(sel_s)), collapse = ", "),
                  N       = paste(sort(unique(N)), collapse = ", "),
                  .groups = "drop"),
      width = Inf)

# N axis per approach: the union over all conditions, so every file has the same columns
N_levels <- lapply(split(dat$N, dat$panel), function(x) sort(unique(x)))

# ---------------- block specs ----------------
specs <- list(
  A = list(title = "A  Outliers: one species",
           def   = "Power = share of selected windows in the top 1% of 10,000",
           xlab  = "Selected windows in the genome (N)",
           null  = 0.01),
  B = list(title = "B  Outliers: parallel across three species",
           def   = "Power = share of parallel genes in all three top-1% lists",
           xlab  = "Genes selected in all three species (N)",
           null  = 0),
  C = list(title = "C  Pathway enrichment: one species",
           def   = sprintf("Power = P(10-window pathway ranks in top %s%%)", gsea_top),
           xlab  = "Selected windows in the pathway (N, of 10)",
           null  = as.numeric(gsea_top) / 100),
  D = list(title = "D  Pathway enrichment: parallel across three species",
           def   = sprintf("Power = P(blood coagulation in top %s%% by mean NES)", gsea_top),
           xlab  = sprintf("Selected blood-coagulation genes (N, of %d)", n_term),
           null  = as.numeric(gsea_top) / 100)
)

theme_heat <- function() {
  theme_bw(base_size = text_size) +
    theme(panel.grid       = element_blank(),
          panel.spacing    = unit(0.5, "lines"),
          strip.background = element_rect(fill = "grey90"),
          strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
          axis.text        = element_text(colour = "#6b6a66", size = text_size * 0.8),
          plot.title       = element_text(face = "bold", size = text_size),
          plot.subtitle    = element_text(colour = "#6b6a66", size = text_size * 0.85),
          legend.position  = "right")
}

# One block of the figure: one approach, rows = Gen, columns = f0
make_block <- function(d, k) {
  sp <- specs[[k]]
  d <- d %>% mutate(N_lab = factor(as.character(N), levels = as.character(N_levels[[k]])))

  p <- ggplot(d, aes(N_lab, s_lab, fill = power)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(Gen_lab ~ f0_lab, drop = FALSE) +
    scale_fill_gradientn(colours = fill_ramp, limits = c(0, 1),
                         breaks = seq(0, 1, 0.25), name = "Power") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(drop = FALSE) +
    labs(title = sp$title,
         subtitle = sprintf("%s. Expected with no selection: %g", sp$def, sp$null),
         x = sp$xlab, y = "Selection coefficient (s)") +
    theme_heat()

  if (!is.na(power_goal)) {
    hi <- d %>% filter(power >= power_goal)
    if (nrow(hi) > 0)
      p <- p + geom_tile(data = hi, fill = NA, colour = "#e69f00", linewidth = 0.6)
  }

  if (show_values) {
    p <- p + geom_text(aes(label = sprintf("%.2f", power),
                           colour = power > 0.55),
                       size = text_size * 0.18, show.legend = FALSE) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "#3d3d3a"), guide = "none")
  }
  p
}

empty_block <- function(k, why) {
  ggplot() + labs(title = specs[[k]]$title, subtitle = why) + theme_void(base_size = text_size) +
    theme(plot.title = element_text(face = "bold", size = text_size),
          plot.subtitle = element_text(colour = "#6b6a66", size = text_size * 0.85))
}

# ---------------- plot: one file per Ne x decline ----------------
combos <- dat %>% distinct(Ne, decline) %>% arrange(ne_num(Ne), desc(decline))

for (i in seq_len(nrow(combos))) {
  ne <- combos$Ne[i]; dec <- combos$decline[i]

  dg <- dat %>% filter(Ne == ne, abs(decline - dec) < 1e-9)

  blocks <- lapply(c("A", "B", "C", "D"), function(k) {
    dk <- dg %>% filter(panel == k)
    if (nrow(dk) == 0) empty_block(k, "No data for this Ne / decline")
    else make_block(dk, k)
  })

  fig <- wrap_plots(blocks, ncol = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("Power across parameter space: Ne = %s, %s",
                      ne, ifelse(dec == 1, "constant size", paste0("decline ", dec))),
      subtitle = paste0("Rows: generation time. Columns: starting frequency f0. ",
                        if (!is.na(power_goal))
                          sprintf("Outlined tiles: power >= %g.", power_goal) else ""),
      theme = theme(plot.title = element_text(face = "bold", size = text_size + 3),
                    plot.subtitle = element_text(colour = "#6b6a66", size = text_size)))

  # width grows with the widest N sweep x number of f0 columns
  n_x <- max(lengths(N_levels))
  w <- max(10, 0.35 * n_x * length(f0_levels) + 3)
  h <- max(10, 0.32 * length(s_levels) * length(gen_levels) * 4 + 5)

  stem <- sprintf("power_heatmap.Ne_%s.decline_%.2f", ne, dec)
  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, paste0(stem, ".", ext)), fig,
           width = w, height = h, dpi = 300, bg = "white", limitsize = FALSE)
  }
  message("Saved: ", file.path(out_dir, paste0(stem, ".png")))
}

message("\nDone. Heatmaps written to: ", out_dir)
