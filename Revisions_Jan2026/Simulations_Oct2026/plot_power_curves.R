#!/usr/bin/env Rscript
# ============================================================
# Power curves from D1_RankPower output
#
# Power = share of permutations in which a 10-window pathway with X selected
# windows ranks in the top 5% (or top 10%) of pathways.
#
# Output, for each threshold (top 5%, top 10%) and each generation time:
#
#   output/power_curves/top05pct/Gen_5/
#       power_curves.main.png             baseline Ne / no decline, one panel per f0
#       power_curves.all_conditions.png   rows: Ne > decline, columns: f0
#       Ne_050k.decline_1.00.png          one file per Ne x decline (sorted Ne, decline)
#       ...
#   output/power_curves/top05pct/Gen_10/ ...
#   output/power_curves/top10pct/Gen_5/  ...
#   output/power_curves/top10pct/Gen_10/ ...
#
# Run from the project directory (where output/ lives):
#   Rscript plot_power_curves.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---------------- settings ----------------
results_file <- "output/fst_rank_by_selected_windows.nnull10000.tsv"
out_root     <- "output/power_curves"

base_Ne      <- "100k"   # baseline for the main figure
base_decline <- 1        # 1 = no decline

# Text size: every label scales from these
size_main  <- 13
size_all   <- 10
size_combo <- 12

metrics <- tibble::tibble(
  column = c("prop_top_5pct", "prop_top_10pct"),
  tag    = c("top05pct",      "top10pct"),
  label  = c("Power (pathway ranks in top 5%)",
             "Power (pathway ranks in top 10%)"),
  null   = c(0.05, 0.10)    # expected value for an all-neutral pathway
)

# One blue ramp, light = weak selection, dark = strong selection
s_ramp <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")

# ---------------- data ----------------
res <- read.delim(results_file)

ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000,
                             as.numeric(x))

res <- res %>%
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

s_cols <- setNames(s_ramp[seq_len(nlevels(res$s))], levels(res$s))
gens   <- sort(unique(res$Gen))

# ---------------- shared plot pieces ----------------
curve_layers <- function(base_size, line_width, y_label, null_level,
                         x_breaks = c(0, 2, 4, 6, 8, 10),
                         y_breaks = seq(0, 1, 0.25)) list(
  geom_hline(yintercept = null_level, linetype = "dashed",
             colour = "#8a8984", linewidth = 0.35),
  geom_line(linewidth = line_width),
  scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE),
  scale_x_continuous(breaks = x_breaks),
  scale_y_continuous(limits = c(0, 1), breaks = y_breaks),
  labs(x = "Selected windows in a 10-window pathway", y = y_label),
  theme_minimal(base_size = base_size),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
    strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
    strip.text.y     = element_text(angle = 0, hjust = 0),
    axis.text        = element_text(colour = "#6b6a66"),
    plot.title       = element_text(face = "bold"),
    legend.position  = "right"
  )
)

# ---------------- Ne / decline sensitivity check ----------------
for (m in seq_len(nrow(metrics))) {
  col <- metrics$column[m]
  spread <- res %>%
    group_by(Gen, sel_s, af, X_selected) %>%
    summarise(range = diff(range(.data[[col]])), .groups = "drop")
  cat(sprintf("%s: spread across Ne x decline - max %.2f, median %.2f\n",
              col, max(spread$range), median(spread$range)))
}

# ---------------- plots ----------------
for (m in seq_len(nrow(metrics))) {

  col   <- metrics$column[m]
  ylab  <- metrics$label[m]
  nullv <- metrics$null[m]

  for (g in gens) {

    out_dir <- file.path(out_root, metrics$tag[m], paste0("Gen_", g))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    dg <- res %>% filter(Gen == g)

    # 1. Main figure: baseline Ne and decline
    d_main <- dg %>% filter(Ne == base_Ne, decline == base_decline)
    if (nrow(d_main) > 0) {
      p_main <- ggplot(d_main, aes(X_selected, .data[[col]], colour = s, group = s)) +
        curve_layers(size_main, 0.8, ylab, nullv) +
        geom_point(size = 1.4) +
        facet_grid(. ~ f0) +
        labs(title = sprintf("Power to detect a selected pathway | %s generations", g),
             subtitle = sprintf("Ne = %s, no decline. Dashed line: expectation for an all-neutral pathway (%g).",
                                base_Ne, nullv))
      ggsave(file.path(out_dir, "power_curves.main.png"), p_main,
             width = 14, height = 4.2, dpi = 300, bg = "white")
    } else {
      warning("No rows for baseline Ne = ", base_Ne, ", decline = ", base_decline, ", Gen = ", g)
    }

    # 2. All conditions for this generation time: rows Ne > decline, columns f0
    n_rows <- n_distinct(dg$Ne_lab) * n_distinct(dg$Dec_lab)
    p_all <- ggplot(dg, aes(X_selected, .data[[col]], colour = s, group = s)) +
      curve_layers(size_all, 0.55, ylab, nullv,
                   x_breaks = c(0, 5, 10), y_breaks = c(0, 0.5, 1)) +
      facet_grid(Ne_lab + Dec_lab ~ f0) +
      labs(title = sprintf("Power across all conditions | %s generations", g),
           subtitle = sprintf("Dashed line: all-neutral expectation (%g)", nullv))
    ggsave(file.path(out_dir, "power_curves.all_conditions.png"), p_all,
           width = 14, height = 1.7 * n_rows + 1.6, dpi = 300, bg = "white",
           limitsize = FALSE)

    # 3. One file per Ne x decline, named so they sort Ne then decline
    combos <- dg %>%
      distinct(Ne, Ne_num, decline, Dec_lab) %>%
      arrange(Ne_num, desc(decline))

    for (i in seq_len(nrow(combos))) {
      cb <- combos[i, ]
      di <- dg %>% filter(Ne == cb$Ne, decline == cb$decline)

      p_i <- ggplot(di, aes(X_selected, .data[[col]], colour = s, group = s)) +
        curve_layers(size_combo, 0.75, ylab, nullv) +
        geom_point(size = 1.2) +
        facet_grid(. ~ f0) +
        labs(title = sprintf("Ne %s | %s generations | %s", cb$Ne, g, cb$Dec_lab))

      fname <- sprintf("Ne_%sk.decline_%.2f.png",
                       formatC(cb$Ne_num / 1000, width = 3, flag = "0"),  # 050, 100, 200 sort correctly
                       cb$decline)
      ggsave(file.path(out_dir, fname), p_i,
             width = 14, height = 4, dpi = 300, bg = "white")
    }

    message("Saved: ", out_dir)
  }
}

message("\nDone. Plots in ", out_root, "/{top05pct,top10pct}/Gen_*/")
