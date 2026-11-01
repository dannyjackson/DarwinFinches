#!/usr/bin/env Rscript
# ============================================================
# Power heatmaps per approach: all Ne values x f0 in one figure
#
# One figure per approach x decline x generation time.
#   approaches  A  outliers: one species           (N = selected windows in genome)
#               B  outliers: three-species overlap (N = parallel genes in genome)
#               C  GSEA: one species               (N = selected windows in 10-window pathway)
#               D  GSEA: three species             (N = selected blood-coagulation genes)
#
#   rows       = effective population size Ne
#   columns    = starting frequency f0
#   tile       = x: number selected (N), y: selection coefficient s, fill: power
#   outline    = stepped orange outlines for high power (default > 0.80, > 0.90, > 0.95)
#
# Companion to power_heatmap_nsel.byNe.R (same inputs, colours, tick and text
# settings); that script puts the four approaches in one figure per Ne, this
# one puts the Ne values in one figure per approach.
#
# Inputs:
#   A  output/fst_outlier_power_nsel.tsv                    outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv             parallel_outlier_overlap.R with
#                                                           PARALLEL_VALUES <- c(1, 2, 5, 10, 20, 30, 50, 100)
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv           plot_gsea_parallel_nsel.R
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/power_heatmap_nsel.byApproach.R
#
# Writes: output/figure_power_heatmap_nsel_byApproach/
#           power_heatmap.<A_outliers_one_species|...>.decline_<d>.Gen_<g>.png / .pdf
#           power_heatmap.byApproach.data.tsv   (values plotted)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---------------- settings ----------------
file_outlier <- "output/fst_outlier_power_nsel.tsv"
file_overlap <- "output/parallel_overlap.permutations.tsv"
file_gsea1   <- "output/fst_rank_by_selected_windows.nnull10000.tsv"
file_gsea3   <- "gsea_empirical_power_nsel/summary_nsel.tsv"
out_dir      <- "output/figure_power_heatmap_nsel_byApproach"

# NULL = every value present in the inputs; or e.g. c("25k", "50k", "75k")
Ne_values      <- NULL   # rows
decline_values <- NULL   # e.g. 1 for constant size only
f0_values      <- NULL   # columns; e.g. c(0.05, 0.2, 0.6)
s_values       <- NULL   # tile rows; e.g. c(0.05, 0.1, 0.2)
approaches     <- c("A", "B", "C", "D")   # which figures to make

gsea_top    <- "5"       # GSEA cutoff for C and D: "5" or "1" (top %)
# Outline levels: each tile gets the outline of the highest level it exceeds.
# Light -> dark orange, all the same thickness. NULL = no outlines. Up to 4 levels.
power_levels  <- c(0.80, 0.90, 0.95)
outline_cols  <- c("#fdc086", "#f28e2b", "#d95f02", "#8c2d04")
outline_width <- 1.2     # line thickness for every level
show_values <- FALSE     # print the power value in each tile (crowded at print size)
show_defs   <- FALSE     # power definition + null under the title
                         # (FALSE for the main text: put these in the caption)

# Tick labels on the discrete axes:
#   "ends" = first and last value only
#   an integer k = every k-th value, always including the last (e.g. 2)
#   "all"  = every value
x_ticks <- "ends"        # N
y_ticks <- "ends"        # s

# Print size. Text sizes are in points AT THIS SIZE, so set the size close to
# how the figure will appear in the paper (e.g. 7.2 in = two-column width).
# At 12 in shrunk to 7.2 in wide, 14 pt -> ~8.5 pt.
fig_width  <- 12         # inches
fig_height <- 7
text_size  <- 14         # base size; titles, strips and ticks scale from this

# light = low power, dark = high power
fill_ramp <- c("#f7f9fc", "#b7d3f6", "#6da7ec", "#2a78d6", "#104281", "#0d366b")

# ---------------- helpers ----------------
ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000,
                             as.numeric(x))

# Stepped outlines: each tile gets the highest power level it exceeds,
# drawn low -> high so the darkest outline sits on top where neighbours overlap
add_outlines <- function(p, d) {
  if (length(power_levels) == 0) return(p)
  lv  <- sort(power_levels)
  if (length(lv) > length(outline_cols)) stop("More power_levels than outline_cols")
  lab <- paste0("> ", formatC(lv, format = "f", digits = 2))
  hi <- d %>%
    mutate(level = findInterval(power, lv, left.open = TRUE)) %>%
    filter(level > 0) %>%
    arrange(level) %>%
    mutate(level_lab = factor(lab[level], levels = lab))
  p +
    geom_tile(data = hi, aes(colour = level_lab), fill = NA,
              linewidth = outline_width, show.legend = TRUE) +
    scale_colour_manual(values = setNames(outline_cols[seq_along(lv)], lab),
                        limits = lab, drop = FALSE, name = "Outline\n(power)") +
    guides(colour = guide_legend(order = 2, override.aes = list(fill = "white")),
           fill = guide_colourbar(order = 1))
}

# Which levels of a discrete axis get a tick label
thin_breaks <- function(lv, mode) {
  n <- length(lv)
  if (n <= 2 || identical(mode, "all")) return(lv)
  idx <- if (identical(mode, "ends")) c(1, n) else unique(c(seq(1, n, by = as.integer(mode)), n))
  lv[idx]
}

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
dat <- bind_rows(A, B, C, D) %>% filter(panel %in% approaches)

if (!is.null(Ne_values))      dat <- dat %>% filter(Ne %in% Ne_values)
if (!is.null(decline_values)) dat <- dat %>% filter(round(decline, 6) %in% round(decline_values, 6))
if (!is.null(f0_values))      dat <- dat %>% filter(round(f0, 6) %in% round(f0_values, 6))
if (!is.null(s_values))       dat <- dat %>% filter(round(sel_s, 6) %in% round(s_values, 6))

if (nrow(dat) == 0) stop("No rows left after filtering; check the settings.")

# Fixed levels so every figure uses the same rows, columns and order
Ne_levels <- unique(dat$Ne[order(ne_num(dat$Ne))])
f0_levels <- sort(unique(dat$f0))
s_levels  <- sort(unique(dat$sel_s))

dat <- dat %>%
  mutate(Ne_lab = factor(paste0("Ne = ", Ne), levels = paste0("Ne = ", Ne_levels)),
         f0_lab = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", f0_levels)),
         s_lab  = factor(as.character(sel_s), levels = as.character(s_levels)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(panel, Ne, decline, Gen, f0, sel_s, N, power) %>%
              arrange(panel, decline, Gen, ne_num(Ne), f0, sel_s, N),
            file.path(out_dir, "power_heatmap.byApproach.data.tsv"),
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

# ---------------- approach specs ----------------
specs <- list(
  A = list(title = "Outliers: one species",
           slug  = "A_outliers_one_species",
           def   = "Power = share of selected windows in the top 1% of 10,000",
           xlab  = "Selected windows (N)",
           null  = 0.01),
  B = list(title = "Outliers: parallel across three species",
           slug  = "B_outliers_three_species",
           def   = "Power = share of parallel genes in all three top-1% lists",
           xlab  = "Parallel genes (N)",
           null  = 0),
  C = list(title = "Pathway enrichment: one species",
           slug  = "C_gsea_one_species",
           def   = sprintf("Power = P(10-window pathway ranks in top %s%%)", gsea_top),
           xlab  = "Selected windows in pathway (N)",
           null  = as.numeric(gsea_top) / 100),
  D = list(title = "Pathway enrichment: parallel across three species",
           slug  = "D_gsea_three_species",
           def   = sprintf("Power = P(blood coagulation in top %s%% by mean NES)", gsea_top),
           xlab  = sprintf("Selected genes in term (N, of %d)", n_term),
           null  = as.numeric(gsea_top) / 100)
)

theme_heat <- function() {
  theme_bw(base_size = text_size) +
    theme(panel.grid        = element_blank(),
          panel.spacing     = unit(0.4, "lines"),
          strip.background  = element_rect(fill = "grey90", colour = NA),
          strip.text        = element_text(face = "bold", colour = "#3d3d3a",
                                           size = text_size * 0.85,
                                           margin = margin(2, 2, 2, 2)),
          axis.text         = element_text(colour = "#3d3d3a", size = text_size * 0.75),
          axis.title        = element_text(size = text_size * 0.85),
          axis.ticks        = element_line(linewidth = 0.4),
          plot.title        = element_text(face = "bold", size = text_size * 1.1),
          plot.subtitle     = element_text(colour = "#6b6a66", size = text_size * 0.75),
          legend.position   = "right",
          legend.title      = element_text(size = text_size * 0.85),
          legend.text       = element_text(size = text_size * 0.75),
          legend.key.height = unit(text_size * 0.12, "lines"))
}

# ---------------- plot: one file per approach x decline x Gen ----------------
combos <- dat %>% distinct(panel, decline, Gen) %>% arrange(panel, desc(decline), Gen)

for (i in seq_len(nrow(combos))) {
  k <- combos$panel[i]; dec <- combos$decline[i]; g <- combos$Gen[i]
  sp <- specs[[k]]

  d <- dat %>%
    filter(panel == k, abs(decline - dec) < 1e-9, Gen == g) %>%
    mutate(N_lab = factor(as.character(N), levels = as.character(N_levels[[k]])))

  missing_ne <- setdiff(Ne_levels, unique(d$Ne))
  if (length(missing_ne) > 0)
    warning(sprintf("%s, decline %g, Gen %g: no data for Ne = %s",
                    sp$title, dec, g, paste(missing_ne, collapse = ", ")))

  subtitle <- paste0(
    sprintf("%s, %s generations. ",
            ifelse(dec == 1, "Constant size", paste0("Decline ", dec)), g),
    if (length(power_levels) > 0) "Orange outlines: high power (see legend)." else "",
    if (show_defs) sprintf("\n%s. Expected with no selection: %g.", sp$def, sp$null) else "")

  p <- ggplot(d, aes(N_lab, s_lab, fill = power)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    facet_grid(Ne_lab ~ f0_lab, drop = FALSE) +
    scale_fill_gradientn(colours = fill_ramp, limits = c(0, 1),
                         breaks = seq(0, 1, 0.25), name = "Power") +
    scale_x_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(N_levels[[k]]), x_ticks)) +
    scale_y_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(s_levels), y_ticks)) +
    labs(title = sp$title, subtitle = subtitle, x = sp$xlab, y = "s") +
    theme_heat()

  p <- add_outlines(p, d)

  if (show_values) {
    p <- p +
      geom_text(data = filter(d, power >  0.55), aes(label = sprintf("%.2f", power)),
                colour = "white",   size = text_size * 0.2) +
      geom_text(data = filter(d, power <= 0.55), aes(label = sprintf("%.2f", power)),
                colour = "#3d3d3a", size = text_size * 0.2)
  }

  stem <- sprintf("power_heatmap.%s.decline_%.2f.Gen_%s", sp$slug, dec, g)
  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, paste0(stem, ".", ext)), p,
           width = fig_width, height = fig_height, dpi = 300, bg = "white")
  }
  message("Saved: ", file.path(out_dir, paste0(stem, ".png")))
}

message("\nDone. Heatmaps written to: ", out_dir)
