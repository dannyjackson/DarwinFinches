#!/usr/bin/env Rscript
# ============================================================
# Composite main-text figure: power heatmaps
#
# One figure per generation time. Two parts stacked, I above II; columns (f0)
# line up between the parts.
#
#   I   Comparison across approaches at one demography (default Ne = 25k, decline = 1)
#         rows    = approach A-D
#         columns = starting frequency f0
#   II  Selected approaches (default D only) at the other Ne values
#       (default 50k, 75k), same decline
#         rows    = one sub-row per Ne within each approach
#         columns = starting frequency f0
#
#   Approaches  A  outliers: one species           (N = selected windows in genome)
#               B  outliers: three-species overlap (N = parallel genes in genome)
#               C  GSEA: one species               (N = selected windows in 10-window pathway)
#               D  GSEA: three species             (N = selected blood-coagulation genes)
#
#   tile    = x: number selected (N), y: selection coefficient s, fill: power
#   outline = orange, one shade per power level (default > 0.80, > 0.90, > 0.95)
#
# No figure title or subtitle: Ne, decline, Gen, power definitions
# and nulls go in the caption (printed to the console at the end as a reminder).
#
# Inputs (same as power_heatmap_nsel.byNe.R):
#   A  output/fst_outlier_power_nsel.tsv                    outlier_power_nsel.R
#   B  output/parallel_overlap.permutations.tsv             parallel_outlier_overlap.R
#   C  output/fst_rank_by_selected_windows.nnull10000.tsv   D1_RankPower.top5.test.R
#   D  gsea_empirical_power_nsel/summary_nsel.tsv           plot_gsea_parallel_nsel.R
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/figure_power_composite.R
#
# Writes: output/figure_power_composite/
#           power_composite.Gen_<g>.png / .pdf
#           power_composite.data.tsv   (values plotted)
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
out_dir      <- "output/figure_power_composite"

decline_main <- 1                       # both parts
Ne_main      <- "25k"                   # part I
Ne_compare   <- c("50k", "75k")         # part II sub-rows (Ne_main is not repeated)
approaches_II <- c("E")                 # approaches shown in part II, e.g. c("C", "D")

# Headings next to the numerals (sentence case, no final period). NULL = numeral only.
# Update the Ne in title_I if you change Ne_main.
title_I  <- "Comparison of approaches (Ne = 25,000)"
title_II <- "Effect of Ne on cross-species pathway enrichment"
f0_values    <- NULL                    # columns in both parts; NULL = every f0 present
Gen_values   <- NULL                    # one figure per Gen; NULL = every Gen present
s_values     <- NULL                    # tile rows; NULL = every s present

gsea_top <- "5"                         # GSEA cutoff for C and D: "5" or "1" (top %)

# Outline levels: each tile gets the outline of the highest level it exceeds.
# Light -> dark orange, all the same thickness. NULL = no outlines. Up to 4 levels.
power_levels  <- c(0.80, 0.90, 0.95)
outline_cols  <- c("#fdc086", "#f28e2b", "#d95f02", "#8c2d04")
outline_width <- 0.45                   # mm, at final size
tile_border   <- 0.1                    # white line between tiles, mm

# Tick labels on the discrete axes: "ends", an integer k (every k-th), or "all"
x_ticks <- "ends"                       # N
y_ticks <- "ends"                       # s

# Final print size (Science / Science Advances). Drawn at 100%, so text sizes
# are the sizes that will appear in the journal.
#   Science widths: 5.7 cm (1 column), 12.1 cm (2 columns), 18.4 cm (3 columns)
#   Text ~7 pt at final size, nothing below 5 pt.
# Height: keep within the printed page (about 22-23 cm); check the current
# Science / Science Advances figure guide before submission.
fig_width_cm  <- 12.1
fig_height_cm <- 21
text_size     <- 7                      # pt: titles; ticks and legend text are 1 pt smaller

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

prep <- function(d, f0_col) {
  d %>%
    mutate(Ne = as.character(Ne), Gen = as.numeric(Gen),
           decline = as.numeric(decline), sel_s = as.numeric(sel_s),
           f0 = as.numeric(.data[[f0_col]])) %>%
    filter(sel_s > 0, abs(decline - decline_main) < 1e-9)
}

# ---------------- data ----------------
A <- read.delim(file_outlier) %>% prep("af") %>%
  transmute(panel = "A", Ne, Gen, sel_s, f0, N = n_selected, power = recall)

B <- read.delim(file_overlap) %>% prep("af") %>%
  filter(n_parallel > 0) %>%
  group_by(Ne, Gen, sel_s, f0, N = n_parallel) %>%
  summarise(power = mean(parallel_in_all_tops), .groups = "drop") %>%
  mutate(panel = "B")

C <- read.delim(file_gsea1) %>% prep("af") %>%
  transmute(panel = "C", Ne, Gen, sel_s, f0, N = X_selected,
            power = .data[[paste0("prop_top_", gsea_top, "pct")]])

D_raw <- read.delim(file_gsea3)
D <- D_raw %>% prep("f0") %>%
  transmute(panel = "D", Ne, Gen, sel_s, f0, N = n_selected,
            power = .data[[paste0("power_cross_top", gsea_top)]])
n_term <- max(D_raw$n_selected)

dat <- bind_rows(A, B, C, D)
if (!is.null(Gen_values)) dat <- dat %>% filter(Gen %in% Gen_values)
if (!is.null(s_values))   dat <- dat %>% filter(round(sel_s, 6) %in% round(s_values, 6))

if (!is.null(f0_values)) dat <- dat %>% filter(round(f0, 6) %in% round(f0_values, 6))

dat_I  <- dat %>% filter(Ne == Ne_main)
dat_II <- dat %>% filter(Ne %in% setdiff(Ne_compare, Ne_main), panel %in% approaches_II)

if (nrow(dat_I) == 0)  stop("Part I: no data for Ne = ", Ne_main, ", decline = ", decline_main)
if (nrow(dat_II) == 0) stop("Part II: no data for Ne = ", paste(Ne_compare, collapse = ", "),
                            ", decline = ", decline_main)

# Shared levels so both parts use the same tiles, columns and colours
s_levels  <- sort(unique(c(dat_I$sel_s, dat_II$sel_s)))
f0_levels <- sort(unique(c(dat_I$f0, dat_II$f0)))
Ne_levels <- unique(dat_II$Ne[order(ne_num(dat_II$Ne))])
N_levels  <- lapply(split(c(dat_I$N, dat_II$N), c(dat_I$panel, dat_II$panel)),
                    function(x) sort(unique(x)))

lab_f0 <- function(d) d %>%
  mutate(col_lab = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", f0_levels)))
dat_I  <- lab_f0(dat_I)
dat_II <- lab_f0(dat_II) %>%
  mutate(Ne_lab = factor(paste0("Ne = ", Ne), levels = paste0("Ne = ", Ne_levels)))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.table(bind_rows(mutate(dat_I, part = "I"), mutate(dat_II, part = "II")) %>%
              select(part, panel, Gen, Ne, f0, sel_s, N, power) %>%
              arrange(Gen, part, panel, ne_num(Ne), f0, sel_s, N),
            file.path(out_dir, "power_composite.data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Coverage check: every approach should be present in both parts for every Gen
cov <- bind_rows(mutate(dat_I, part = "I"), mutate(dat_II, part = "II")) %>%
  distinct(part, Gen, panel, Ne, f0)
cat("Coverage per part x Gen x approach:\n")
print(cov %>% group_by(part, Gen, panel) %>%
        summarise(Ne = paste(unique(Ne[order(ne_num(Ne))]), collapse = ", "),
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

small_size <- max(5, text_size - 1)     # ticks, strips, legend text (>= 5 pt)

theme_heat <- function() {
  theme_bw(base_size = text_size, base_line_size = 0.25, base_rect_size = 0.25) +
    theme(panel.grid        = element_blank(),
          panel.spacing     = unit(1.2, "mm"),
          strip.background  = element_rect(fill = "grey90", colour = NA),
          strip.text        = element_text(face = "bold", colour = "#3d3d3a",
                                           size = small_size,
                                           margin = margin(1, 1, 1, 1)),
          axis.text         = element_text(colour = "#3d3d3a", size = small_size),
          axis.title        = element_text(size = text_size),
          axis.title.y      = element_text(margin = margin(0, 1, 0, 0)),
          axis.title.x      = element_text(margin = margin(1, 0, 0, 0)),
          axis.ticks        = element_line(linewidth = 0.25),
          axis.ticks.length = unit(0.6, "mm"),
          plot.title        = element_text(face = "bold", size = text_size,
                                           margin = margin(0, 0, 1.5, 0)),
          plot.title.position = "plot",
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

# One row of one part: one approach, columns = f0; part II adds one sub-row per Ne
make_row <- function(d, k, top, part) {
  sp <- specs[[k]]
  if (nrow(d) == 0) return(plot_spacer())
  d <- d %>%
    mutate(N_lab = factor(as.character(N), levels = as.character(N_levels[[k]])),
           s_lab = factor(as.character(sel_s), levels = as.character(s_levels)))

  p <- ggplot(d, aes(N_lab, s_lab, fill = power)) +
    geom_tile(colour = "white", linewidth = tile_border) +
    (if (part == "I") facet_wrap(~ col_lab, nrow = 1, drop = FALSE)
     else facet_grid(Ne_lab ~ col_lab, drop = FALSE)) +
    scale_fill_gradientn(colours = fill_ramp, limits = c(0, 1),
                         breaks = seq(0, 1, 0.25), name = "Power") +
    scale_x_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(N_levels[[k]]), x_ticks)) +
    scale_y_discrete(drop = FALSE,
                     breaks = thin_breaks(as.character(s_levels), y_ticks)) +
    labs(title = sp$title, x = sp$xlab, y = "s") +
    theme_heat()

  # f0 column labels on the top row of each part only; Ne row labels (II) on every row
  if (!top) p <- p + theme(strip.text.x = element_blank(),
                           strip.background.x = element_blank())

  add_outlines(p, d)
}

# Part header: bold numeral, then a bold heading after a small gap
numeral <- function(txt, title = NULL) {
  num <- grid::textGrob(txt, x = grid::unit(0, "npc"), hjust = 0,
                        gp = grid::gpar(fontsize = text_size + 2, fontface = "bold"))
  if (is.null(title) || !nzchar(title)) return(wrap_elements(full = num))
  ttl <- grid::textGrob(title,
                        x = grid::unit(0, "npc") + grid::grobWidth(num) + grid::unit(2, "mm"),
                        hjust = 0,
                        gp = grid::gpar(fontsize = text_size + 1, fontface = "bold"))
  wrap_elements(full = grid::grobTree(num, ttl))
}

# ---------------- plot: one figure per Gen ----------------
for (g in sort(unique(c(dat_I$Gen, dat_II$Gen)))) {
  dI  <- dat_I  %>% filter(Gen == g)
  dII <- dat_II %>% filter(Gen == g)

  ks    <- c("A", "B", "C", "D")
  ks_II <- intersect(ks, approaches_II)
  cells <- c(list(numeral("I", title_I)),
             lapply(ks, function(k) make_row(filter(dI, panel == k), k, k == ks[1], "I")),
             list(numeral("II", title_II)),
             lapply(ks_II, function(k) make_row(filter(dII, panel == k), k,
                                                k == ks_II[1], "II")))

  n_ne <- max(1, length(Ne_levels))
  fig <- wrap_plots(cells, ncol = 1,
                    heights = c(0.08, rep(1, 4), 0.08, rep(n_ne, length(ks_II)))) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.margin = margin(0, 4, 0, 4))

  stem <- sprintf("power_composite.Gen_%s", g)
  for (ext in c("png", "pdf")) {
    ggsave(file.path(out_dir, paste0(stem, ".", ext)), fig,
           width = fig_width_cm, height = fig_height_cm, units = "cm",
           dpi = 600, bg = "white")
  }
  message("Saved: ", file.path(out_dir, paste0(stem, ".png")))
}

# ---------------- caption reminder ----------------
cat("\nFor the caption:\n")
cat(sprintf("  I: Ne = %s, approaches A-D. II: approach %s at Ne = %s. Both %s; columns = f0.\n",
            Ne_main, paste(approaches_II, collapse = ", "), paste(Ne_levels, collapse = ", "),
            ifelse(decline_main == 1, "constant size", paste0("decline ", decline_main))))
cat("  A: power = share of selected windows in the top 1% of 10,000 (null 0.01)\n")
cat("  B: power = share of parallel genes in all three top-1% lists (null 0)\n")
cat(sprintf("  C: power = P(10-window pathway ranks in top %s%%) (null %g)\n",
            gsea_top, as.numeric(gsea_top) / 100))
cat(sprintf("  D: power = P(blood coagulation in top %s%% by mean NES) (null %g)\n",
            gsea_top, as.numeric(gsea_top) / 100))
if (length(power_levels) > 0)
  cat(sprintf("  Orange outlines: power > %s (light to dark)\n",
              paste(formatC(sort(power_levels), format = "f", digits = 2), collapse = ", ")))
