#!/usr/bin/env Rscript
# ============================================================
# Outlier-window power analysis + power curves
#
# Replaces E_Simulations/Power_analysis_two_time_points/plot.power_analysis.r
# for the current simulations (Simulations_Oct2026 output layout).
#
# Design, per selected condition (Ne x Gen x decline x s x f0):
#   - draw N_NEUTRAL (9,900) windows from the matching neutral file
#     (same Ne, Gen, decline; d*_s0.0_af0.5.summary_stats.tsv)
#   - draw N_SELECTED (100) windows from the selected file
#   - rank all 10,000 windows by FST
#   - power = proportion of the top TOP_N (100) windows that are selected
#   - repeat N_PERMUTATIONS (100) times and summarise
#
# Because N_SELECTED == TOP_N, this proportion is also the fraction of
# selected windows that reach the top 1% (precision = recall here).
# With no power, the expected value is N_SELECTED / N_TOTAL = 0.01.
#
# Run from the project directory, e.g.
#   cd /xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
#   Rscript ${SCRIPTDIR}/outlier_power_curves.R
#
# Writes:
#   output/fst_outlier_power.nnull9900.nsel100.tsv        (summary table)
#   output/fst_outlier_power.permutations.tsv             (every permutation)
#   output/power_curves_outlier/top1pct/Gen_<g>/
#       power_curves.main.png             Ne = 100k, no decline
#       power_curves.all_conditions.png   Ne x decline grid
#       Ne_<NNN>k.decline_<d>.png         one figure per Ne x decline
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

set.seed(123)

# ============================================================
# Settings
# ============================================================

root_dir <- "."

N_NEUTRAL      <- 9900
N_SELECTED     <- 100
N_TOTAL        <- N_NEUTRAL + N_SELECTED
TOP_N          <- 100     # top 1% of 10,000 windows
N_PERMUTATIONS <- 100

NULL_POWER <- N_SELECTED / N_TOTAL   # expected value with no power (0.01)

out_dir  <- "output"
out_root <- file.path(out_dir, "power_curves_outlier", "top1pct")

base_Ne      <- "100k"   # baseline for the main figure
base_decline <- 1        # 1 = no decline

# Shade the 2.5-97.5% range of power across permutations (FALSE = lines only,
# matching the GSEA curves)
SHOW_RANGE <- FALSE

size_main  <- 13
size_all   <- 10
size_combo <- 12

# Same blue ramp as plot_power_curves.R (one colour per s value)
s_ramp <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")

# ============================================================
# Power for one selected condition
# ============================================================

outlier_power <- function(null_fst, sel_fst, n_permutations) {
  is_sel <- c(rep(FALSE, N_NEUTRAL), rep(TRUE, N_SELECTED))

  vapply(seq_len(n_permutations), function(perm) {
    fst <- c(sample(null_fst, N_NEUTRAL), sample(sel_fst, N_SELECTED))
    # Shuffle before ordering so ties are broken at random
    shuffle <- sample.int(N_TOTAL)
    top <- shuffle[order(-fst[shuffle])][seq_len(TOP_N)]
    mean(is_sel[top])
  }, numeric(1))
}

read_fst <- function(file) {
  d <- read.delim(file)
  d$fst_hudson[is.finite(d$fst_hudson)]
}

# ============================================================
# Find files (same layout as D1_RankPower.top5.test.R)
# ============================================================

all_files <- list.files(
  root_dir,
  pattern = "summary_stats\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

all_files <- all_files[
  str_detect(all_files, "^\\./Ne_[^/]+/Gen_[0-9]+/[^/]+$")
]

selected_files <- all_files[!str_detect(basename(all_files), "_s0\\.0_")]

if (length(selected_files) == 0) {
  stop("No selected *.summary_stats.tsv files found under ", root_dir,
       " (run this from the gsea_power_analysis directory)")
}

# Cache neutral FST so each neutral file is read once
neutral_cache <- list()

get_neutral <- function(selected_file, decline) {
  selected_dir <- dirname(selected_file)

  neutral_candidates <- list.files(
    selected_dir,
    pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$",
    full.names = TRUE
  )
  if (length(neutral_candidates) == 0) return(NULL)

  candidate_declines <- as.numeric(
    str_match(basename(neutral_candidates), "^d([0-9.]+)_s")[, 2]
  )
  match_index <- which(candidate_declines == as.numeric(decline))
  if (length(match_index) == 0) return(NULL)

  null_file <- neutral_candidates[match_index[1]]
  if (is.null(neutral_cache[[null_file]])) {
    neutral_cache[[null_file]] <<- read_fst(null_file)
  }
  list(file = null_file, fst = neutral_cache[[null_file]])
}

# ============================================================
# Run every selected condition
# ============================================================

perm_results <- map_dfr(selected_files, function(selected_file) {

  filename <- basename(selected_file)
  decline  <- str_match(filename, "^d([0-9.]+)_s")[, 2]
  sel_s    <- str_match(filename, "_s([0-9.]+)_af")[, 2]
  af       <- str_match(filename, "_af([0-9.]+)\\.summary")[, 2]
  Ne       <- str_match(selected_file, "Ne_([^/]+)")[, 2]
  Gen      <- str_match(selected_file, "Gen_([0-9]+)")[, 2]

  neutral <- get_neutral(selected_file, decline)
  if (is.null(neutral)) {
    warning("No matching neutral file for: ", selected_file)
    return(NULL)
  }
  if (length(neutral$fst) < N_NEUTRAL) {
    warning("Not enough usable neutral windows: ", neutral$file,
            " n = ", length(neutral$fst))
    return(NULL)
  }

  sel_fst <- read_fst(selected_file)
  if (length(sel_fst) < N_SELECTED) {
    warning("Not enough usable selected windows: ", selected_file,
            " n = ", length(sel_fst))
    return(NULL)
  }

  cat("Processing: Ne =", Ne, "Gen =", Gen, "d =", decline,
      "s =", sel_s, "af =", af,
      "| neutral n =", length(neutral$fst), "selected n =", length(sel_fst), "\n")

  power <- outlier_power(neutral$fst, sel_fst, N_PERMUTATIONS)

  tibble(
    Ne = Ne,
    Gen = as.numeric(Gen),
    decline = as.numeric(decline),
    sel_s = as.numeric(sel_s),
    af = as.numeric(af),
    n_neutral_pool = length(neutral$fst),
    n_selected_pool = length(sel_fst),
    permutation = seq_along(power),
    power = power
  )
})

if (nrow(perm_results) == 0) stop("No conditions could be analysed.")

results <- perm_results %>%
  group_by(Ne, Gen, decline, sel_s, af, n_neutral_pool, n_selected_pool) %>%
  summarise(
    n_permutations = n(),
    mean_power   = mean(power),
    median_power = median(power),
    sd_power     = sd(power),
    q025_power   = quantile(power, 0.025),
    q975_power   = quantile(power, 0.975),
    .groups = "drop"
  ) %>%
  mutate(null_power = NULL_POWER)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results_file <- file.path(
  out_dir,
  paste0("fst_outlier_power.nnull", N_NEUTRAL, ".nsel", N_SELECTED, ".tsv")
)
write.table(results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(perm_results, file.path(out_dir, "fst_outlier_power.permutations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", results_file)

# ============================================================
# Power curves (format follows plot_power_curves.R)
#   x = selection coefficient, colour = s, facets = f0
# ============================================================

ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000,
                             as.numeric(x))

res <- results %>%
  mutate(
    Ne      = as.character(Ne),
    Ne_num  = ne_num(Ne),
    s       = factor(sel_s, levels = sort(unique(sel_s))),
    f0      = factor(paste0("f0 = ", af),
                     levels = paste0("f0 = ", sort(unique(af)))),
    Ne_lab  = factor(paste0("Ne ", Ne),
                     levels = paste0("Ne ", unique(Ne[order(Ne_num)]))),
    Dec_lab = factor(ifelse(decline == 1, "no decline", paste0("decline ", decline)),
                     levels = c("no decline",
                                paste0("decline ",
                                       sort(unique(decline[decline != 1]), decreasing = TRUE))))
  )

s_cols <- setNames(rep_len(s_ramp, nlevels(res$s)), levels(res$s))
s_vals <- sort(unique(res$sel_s))
gens   <- sort(unique(res$Gen))

ylab <- "Power (selected windows in top 1%)"

curve_layers <- function(base_size, line_width, y_breaks = seq(0, 1, 0.25),
                         x_labels = TRUE) list(
  if (SHOW_RANGE)
    geom_ribbon(aes(ymin = q025_power, ymax = q975_power),
                fill = "#2a78d6", alpha = 0.15, colour = NA),
  geom_hline(yintercept = NULL_POWER, linetype = "dashed",
             colour = "#8a8984", linewidth = 0.35),
  geom_line(linewidth = line_width, colour = "#3d3d3a"),
  scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE),
  scale_x_log10(breaks = s_vals,
                labels = if (x_labels) as.character(s_vals) else NULL),
  scale_y_continuous(limits = c(0, 1), breaks = y_breaks),
  labs(x = "Selection coefficient (s, log scale)", y = ylab),
  theme_minimal(base_size = base_size),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
    strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
    strip.text.y     = element_text(angle = 0, hjust = 0),
    axis.text        = element_text(colour = "#6b6a66"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    legend.position  = "right"
  )
)

# s is on the x axis, so each panel has one curve; points keep the GSEA
# colour key for s.
for (g in gens) {
  out_dir_g <- file.path(out_root, paste0("Gen_", g))
  dir.create(out_dir_g, recursive = TRUE, showWarnings = FALSE)

  dg <- res %>% filter(Gen == g)

  # ---- Main figure: baseline Ne, no decline ----
  d_main <- dg %>% filter(Ne == base_Ne, decline == base_decline)
  if (nrow(d_main) > 0) {
    p_main <- ggplot(d_main, aes(sel_s, mean_power, group = 1)) +
      curve_layers(size_main, 0.8) +
      geom_point(aes(colour = s), size = 2) +
      facet_grid(. ~ f0) +
      labs(title = sprintf("Power to detect outlier windows | %s generations", g),
           subtitle = sprintf(paste0("Ne = %s, no decline. %d selected among %s windows; ",
                                     "dashed line: expectation with no power (%g)."),
                              base_Ne, N_SELECTED, format(N_TOTAL, big.mark = ","),
                              NULL_POWER))
    ggsave(file.path(out_dir_g, "power_curves.main.png"), p_main,
           width = 14, height = 4.2, dpi = 300, bg = "white")
  } else {
    warning("No rows for baseline Ne = ", base_Ne, ", decline = ", base_decline, ", Gen = ", g)
  }

  # ---- All conditions: Ne x decline rows, f0 columns ----
  n_rows <- n_distinct(dg$Ne_lab) * n_distinct(dg$Dec_lab)
  p_all <- ggplot(dg, aes(sel_s, mean_power, group = 1)) +
    curve_layers(size_all, 0.55, y_breaks = c(0, 0.5, 1)) +
    geom_point(aes(colour = s), size = 1.2) +
    facet_grid(Ne_lab + Dec_lab ~ f0) +
    labs(title = sprintf("Outlier power across all conditions | %s generations", g),
         subtitle = sprintf("Dashed line: expectation with no power (%g)", NULL_POWER))
  ggsave(file.path(out_dir_g, "power_curves.all_conditions.png"), p_all,
         width = 14, height = 1.7 * n_rows + 1.6, dpi = 300, bg = "white",
         limitsize = FALSE)

  # ---- One figure per Ne x decline ----
  combos <- dg %>%
    distinct(Ne, Ne_num, decline, Dec_lab) %>%
    arrange(Ne_num, desc(decline))

  for (i in seq_len(nrow(combos))) {
    cb <- combos[i, ]
    di <- dg %>% filter(Ne == cb$Ne, decline == cb$decline)

    p_i <- ggplot(di, aes(sel_s, mean_power, group = 1)) +
      curve_layers(size_combo, 0.75) +
      geom_point(aes(colour = s), size = 1.8) +
      facet_grid(. ~ f0) +
      labs(title = sprintf("Outlier power | Ne %s | %s generations | %s",
                           cb$Ne, g, cb$Dec_lab))

    fname <- sprintf("Ne_%sk.decline_%.2f.png",
                     formatC(cb$Ne_num / 1000, width = 3, flag = "0"),
                     cb$decline)
    ggsave(file.path(out_dir_g, fname), p_i,
           width = 14, height = 4, dpi = 300, bg = "white")
  }

  message("Saved: ", out_dir_g)
}

message("\nDone. Table: ", results_file, "\nPlots: ", out_root)