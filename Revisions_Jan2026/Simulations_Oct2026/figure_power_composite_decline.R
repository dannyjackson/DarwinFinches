#!/usr/bin/env Rscript
# ============================================================
# Effect of population decline on power, paired by Ne
#
# One figure per generation time. For each approach shown (default D only),
# a grid of heatmaps with
#   rows    = Ne, with constant size and decline as adjacent sub-rows
#             (Ne 25k: constant, decline; Ne 50k: constant, decline; ...)
#   columns = starting frequency f0
#   tile    = x: number selected (N), y: selection coefficient s, fill: power
#   outline = orange, one shade per power level (default > 0.80, > 0.90, > 0.95)
#
#   Approaches  A  outliers: one species           (N = selected windows in genome)
#               B  outliers: three-species overlap (N = parallel genes in genome)
#               C  GSEA: one species               (N = selected windows in 10-window pathway)
#               D  GSEA: three species             (N = selected blood-coagulation genes)
#
# Inputs (same as figure_power_composite.R):
#   A  output/fst_outlier_power_nsel.tsv                    outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv             parallel_outlier_overlap.R
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv           plot_gsea_parallel_nsel.R
#
# Every Ne x decline combination requested must be present in the inputs for
# the approaches shown; the script prints a coverage table and stops otherwise.
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/figure_power_composite_decline.R
#
# Writes: output/figure_power_composite_decline/
#           power_decline.Gen_<g>.png / .pdf
#           power_decline.data.tsv   (values plotted)
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
out_dir      <- "output/figure_power_composite_decline"

Ne_values      <- c("25k", "50k", "75k")   # row blocks, in this order
decline_values <- c(1, 0.9)                # sub-rows within each Ne, in this order
approaches     <- c("D")                   # one stacked block per approach, e.g. c("C", "D")

# Sub-row labels, keyed by decline value
decline_labels <- c("1"   = "constant",
                    "0.9" = "decline 10%/gen",
                    "0.8" = "decline 20%/gen")

N_values_B <- c(1, 2, 5, 10, 20, 50, 100)   # parallel genes shown in panel B (drops 30)

f0_values  <- NULL                      # columns; NULL = every f0 present
Gen_values <- NULL                      # one figure per Gen; NULL = every Gen present
s_values   <- NULL                      # tile rows; NULL = every s present

gsea_top <- "5"                         # GSEA cutoff for C and D: "5" or "1" (top %)

# Outline levels: each tile gets the outline of the highest level it exceeds.
power_levels  <- c(0.80, 0.90, 0.95)
outline_cols  <- c("#fdc086", "#f28e2b", "#d95f02", "#8c2d04")
outline_width <- 0.45                   # mm, at final size
tile_border   <- 0.1                    # white line between tiles, mm

# Tick labels on the discrete axes: "ends", an integer k (every k-th), or "all"
x_ticks      <- "all"                   # N
y_ticks      <- "all"                   # s
dodge_x_over <- 8                       # stagger x labels on two rows when more than this many

# Final print size (Science / Science Advances). Text sizes are final sizes.
fig_width_cm  <- 12.1
fig_height_cm <- 4 + 1.7 * length(Ne_values) * length(decline_values) * length(approaches)
text_size     <- 7                      # pt

# light = low power, dark = high power
fill_ramp <- c("#f7f9fc", "#b7d3f6", "#6da7ec", "#2a78d6", "#104281", "#0d366b")

# ---------------- helpers ----------------
ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000,
                             as.numeric(x))

thin_breaks <- function(lv, mode) {
  n <- length(lv)
  if (n <= 2 || identical(mode, "all")) return(lv)
  idx <- if (identical(mode, "ends")) c(1, n) else unique(c(seq(1, n, by = as.integer(mode)), n))
  lv[idx]
}

dec_lab <- function(x) {
  key <- as.character(x)
  out <- decline_labels[key]
  out[is.na(out)] <- paste0("decline ", key[is.na(out)])
  unname(out)
}

same_dec <- function(x, y) abs(x - y) < 1e-9
in_decs  <- function(x, ys) vapply(x, function(v) any(same_dec(v, ys)), logical(1))

prep <- function(d, f0_col) {
  d %>%
    mutate(Ne = as.character(Ne), Gen = as.numeric(Gen),
           decline = as.numeric(decline), sel_s = as.numeric(sel_s),
           f0 = as.numeric(.data[[f0_col]])) %>%
    filter(sel_s > 0, Ne %in% Ne_values, in_decs(decline, decline_values))
}

# ---------------- data ----------------
A <- read.delim(file_outlier) %>% prep("af") %>%
  transmute(panel = "A", Ne, Gen, decline, sel_s, f0, N = n_selected, power = recall)

B <- read.delim(file_overlap) %>% prep("af") %>%
  filter(n_parallel %in% N_values_B) %>%
  group_by(Ne, Gen, decline, sel_s, f0, N = n_parallel) %>%
  summarise(power = mean(parallel_in_all_tops), .groups = "drop") %>%
  mutate(panel = "B")

C <- read.delim(file_gsea1) %>% prep("af") %>%
  transmute(panel = "C", Ne, Gen, decline, sel_s, f0, N = X_selected,
            power = .data[[paste0("prop_top_", gsea_top, "pct")]])

D_raw <- read.delim(file_gsea3)
D <- D_raw %>% prep("f0") %>%
  transmute(panel = "D", Ne, Gen, decline, sel_s, f0, N = n_selected,
            power = .data[[paste0("power_cross_top", gsea_top)]])
n_term <- max(D_raw$n_selected)

dat <- bind_rows(A, B, C, D) %>% filter(panel %in% approaches)
if (!is.null(Gen_values)) dat <- dat %>% filter(Gen %in% Gen_values)
if (!is.null(s_values))   dat <- dat %>% filter(round(sel_s, 6) %in% round(s_values, 6))
if (!is.null(f0_values))  dat <- dat %>% filter(round(f0, 6) %in% round(f0_values, 6))

if (nrow(dat) == 0) stop("No data for the requested Ne / decline / approaches")

# Every Ne x decline must be present for every approach shown
want <- expand.grid(panel = approaches, Ne = Ne_values, decline = decline_values,
                    stringsAsFactors = FALSE)
have <- dat %>% distinct(panel, Ne, decline)
present <- vapply(seq_len(nrow(want)), function(i)
  any(have$panel == want$panel[i] & have$Ne == want$Ne[i] &
        same_dec(have$decline, want$decline[i])), logical(1))
miss <- want[!present, ]
if (nrow(miss) > 0) {
  print(miss)
  stop("Missing conditions (above). Run the simulations / GSEA jobs for them first.")
}

# Row labels: Ne block, then decline sub-row, in the requested order
row_levels <- as.vector(t(outer(Ne_values, decline_values,
                                function(ne, dc) paste0("Ne = ", ne, "\n", dec_lab(dc)))))

s_levels  <- sort(unique(dat$sel_s))
f0_levels <- sort(unique(dat$f0))
N_levels  <- lapply(split(dat$N, dat$panel), function(x) sort(unique(x)))

dat <- dat %>%
  mutate(col_lab = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", f0_levels)),
         row_lab = factor(paste0("Ne = ", Ne, "\n", dec_lab(decline)), levels = row_levels),
         s_lab   = factor(as.character(sel_s), levels = as.character(s_levels)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(dat %>% select(panel, Gen, Ne, decline, f0, sel_s, N, power) %>%
              arrange(Gen, panel, ne_num(Ne), desc(decline), f0, sel_s, N),
            file.path(out_dir, "power_decline.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Coverage per Gen x approach:\n")
print(dat %>% distinct(Gen, panel, Ne, decline, f0) %>%
        group_by(Gen, panel) %>%
        summarise(Ne = paste(unique(Ne[order(ne_num(Ne))]), collapse = ", "),
                  decline = paste(sort(unique(decline), decreasing = TRUE), collapse = ", "),
                  f0 = paste(sort(unique(f0)), collapse = ", "),
                  .groups = "drop"), n = Inf, width = Inf)

# ---------------- approach labels ----------------
specs <- list(
  A = list(title = "A  Outliers, one species",           xlab = "Selected windows (N)"),
  B = list(title = "B  Outliers, three species",         xlab = "Parallel genes (N)"),
  C = list(title = "C  Pathway enrichment, one species",  xlab = "Selected windows in pathway (N)"),
  D = list(title = "D  Pathway enrichment, three species",
           xlab = sprintf("Selected genes in term (N, of %d)", n_term))
)
if (length(approaches) == 1)   # single block: drop the letter
  specs[[approaches]]$title <- sub("^[A-D]  ", "", specs[[approaches]]$title)

small_size <- max(5, text_size - 1)

theme_heat <- function() {
  theme_bw(base_size = text_size, base_line_size = 0.25, base_rect_size = 0.25) +
    theme(panel.grid        = element_blank(),
          panel.spacing.x   = unit(1.2, "mm"),
          panel.spacing.y   = unit(0.8, "mm"),
          strip.background  = element_rect(fill = "grey90", colour = NA),
          strip.text        = element_text(face = "bold", colour = "#3d3d3a",
                                           size = small_size,
                                           margin = margin(1, 1, 1, 1)),
          strip.text.y      = element_text(angle = 0, hjust = 0,
                                           margin = margin(1, 2, 1, 2)),
          axis.text         = element_text(colour = "#3d3d3a", size = small_size),
          axis.title        = element_text(size = text_size),
          axis.title.y      = element_text(margin = margin(0, 1, 0, 0)),
          axis.title.x      = element_text(margin = margin(1, 0, 0, 0)),
          axis.ticks        = element_line(linewidth = 0.25),
          axis.ticks.length = unit(0.6, "mm"),
          plot.title        = element_text(face = "bold", size = text_size,
                                           margin = margin(0, 0, 1.5, 0)),
          plot.title.position = "plot",
          legend.position   = "bottom",
          legend.box        = "horizontal",
          legend.margin     = margin(0, 4, 0, 4),
          legend.title      = element_text(size = text_size),
          legend.text       = element_text(size = small_size),
          legend.key.size   = unit(3, "mm"),
          plot.margin       = margin(1.5, 2, 1.5, 1.5))
}

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
    guides(colour = guide_legend(order = 2, nrow = 1, title.position = "left",
                                 override.aes = list(fill = "white")),
           fill = guide_colourbar(order = 1, title.position = "left",
                                  barwidth = unit(30, "mm"), barheight = unit(2.5, "mm")))
}

# One block: one approach; rows = Ne x decline, columns = f0
make_block <- function(d, k) {
  sp <- specs[[k]]
  d <- d %>% mutate(N_lab = factor(as.character(N), levels = as.character(N_levels[[k]])))
  x_guide <- if (length(N_levels[[k]]) > dodge_x_over) guide_axis(n.dodge = 2) else guide_axis()

  p <- ggplot(d, aes(N_lab, s_lab, fill = power)) +
    geom_tile(colour = "white", linewidth = tile_border) +
    facet_grid(row_lab ~ col_lab, drop = FALSE) +
    scale_fill_gradientn(colours = fill_ramp, limits = c(0, 1),
                         breaks = seq(0, 1, 0.25), name = "Power") +
    scale_x_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(N_levels[[k]]), x_ticks),
                     guide  = x_guide) +
    scale_y_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(s_levels), y_ticks)) +
    labs(title = sp$title, x = sp$xlab, y = "s") +
    theme_heat()

  add_outlines(p, d)
}

# ---------------- plot: one figure per Gen ----------------
for (g in sort(unique(dat$Gen))) {
  dg <- dat %>% filter(Gen == g)
  blocks <- lapply(approaches, function(k) make_block(filter(dg, panel == k), k))

  fig <- if (length(blocks) == 1) blocks[[1]] else
    wrap_plots(blocks, ncol = 1) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  stem <- sprintf("power_decline.Gen_%s", g)
  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, paste0(stem, ".", ext)), fig,
           width = fig_width_cm, height = fig_height_cm, units = "cm",
           dpi = 600, bg = "white")
  }
  message("Saved: ", file.path(out_dir, paste0(stem, ".png")))
}

# ---------------- caption reminder ----------------
cat("\nFor the caption:\n")
cat(sprintf("  Approach %s. Rows: Ne = %s, each at %s. Columns = f0. One figure per Gen.\n",
            paste(approaches, collapse = ", "), paste(Ne_values, collapse = ", "),
            paste(dec_lab(decline_values), collapse = " and ")))
cat("  A: power = share of selected windows in the top 1% of 10,000 (null 0.01)\n")
cat("  B: power = share of parallel genes in all three top-1% lists (null ~0)\n")
cat(sprintf("  C: power = P(10-window pathway ranks in top %s%%) (null %g)\n",
            gsea_top, as.numeric(gsea_top) / 100))
cat(sprintf("  D: power = P(blood coagulation in top %s%% by mean NES) (null %g)\n",
            gsea_top, as.numeric(gsea_top) / 100))
if (length(power_levels) > 0)
  cat(sprintf("  Orange outlines: power > %s (light to dark)\n",
              paste(formatC(sort(power_levels), format = "f", digits = 2), collapse = ", ")))
