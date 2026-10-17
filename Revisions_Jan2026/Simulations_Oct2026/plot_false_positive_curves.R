#!/usr/bin/env Rscript
# ============================================================
# False positive curves for the three-species overlap test
#
# Reads output/parallel_overlap.permutations.tsv (from parallel_outlier_overlap.R)
# and computes the false-positive summaries itself, so it works with output
# from both versions of that script
# and plots, against selection strength:
#   - P(at least one false positive) in the overlap
#   - false discovery rate (false positives / all genes in the overlap)
# for the 3-way overlap (top list in all three species) and the 2+ overlap
# (top list in at least two species).
#
# A false positive = a gene in the overlap that is NOT selected in all three
# species (neutral everywhere, or selected in one species only).
#
# Dashed lines: what you'd expect by chance alone with no selected genes.
#
# Run from the project directory:
#   Rscript ${SCRIPTDIR}/plot_false_positive_curves.R
#
# Writes output/power_curves_parallel/par<NN>/false_positives/Gen_<g>/
#   fp_curves.main.png              baseline Ne / no decline
#   fp_curves.all_conditions.png    rows: Ne > decline, columns: f0 (P(>=1 FP) only)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ---------------- settings (match parallel_outlier_overlap.R) ----------------
perm_file    <- "output/parallel_overlap.permutations.tsv"   # works with either script version
out_root     <- "output/power_curves_parallel"

N_GENES <- 10000
TOP_N   <- 100

base_Ne      <- "100k"
base_decline <- 1

size_main <- 13
size_all  <- 10

# Two series, fixed colours: 3-way = blue, 2+ = orange
overlap_cols <- c("All three species" = "#2a78d6",
                  "At least two species" = "#eb6834")

# ---------------- chance expectation with no selection ----------------
p  <- TOP_N / N_GENES
lambda3 <- N_GENES * p^3                          # expected chance genes in 3-way overlap
lambda2 <- N_GENES * (3 * p^2 * (1 - p) + p^3)    # ... in 2+ overlap
chance <- tibble(
  overlap = names(overlap_cols),
  metric  = "P(at least one false positive)",
  value   = c(1 - exp(-lambda3), 1 - exp(-lambda2))
)
cat(sprintf("Chance with no selection: P(>=1 FP) = %.3f (3-way), %.3f (2+)\n",
            chance$value[1], chance$value[2]))

# ---------------- data ----------------
ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

pooled <- function(num, den) ifelse(sum(den) > 0, sum(num) / sum(den), NA_real_)

perm <- read.delim(perm_file)
need <- c("overlap3_n", "overlap3_parallel", "overlap2_n", "overlap2_parallel")
miss <- setdiff(need, names(perm))
if (length(miss)) stop("Missing columns in ", perm_file, ": ", paste(miss, collapse = ", "))

res <- perm %>%
  mutate(overlap3_false = overlap3_n - overlap3_parallel,
         overlap2_false = overlap2_n - overlap2_parallel) %>%
  group_by(Ne, Gen, decline, sel_s, af, n_parallel) %>%
  summarise(
    n_permutations  = n(),
    prob_any_false3 = mean(overlap3_false > 0),
    prob_any_false2 = mean(overlap2_false > 0),
    fdr3            = pooled(overlap3_false, overlap3_n),
    fdr2            = pooled(overlap2_false, overlap2_n),
    .groups = "drop"
  ) %>%
  mutate(
    Ne      = as.character(Ne),
    Ne_num  = ne_num(Ne),
    f0      = factor(paste0("f0 = ", af), levels = paste0("f0 = ", sort(unique(af)))),
    Ne_lab  = factor(paste0("Ne ", Ne), levels = paste0("Ne ", unique(Ne[order(Ne_num)]))),
    Dec_lab = factor(ifelse(decline == 1, "no decline", paste0("decline ", decline)),
                     levels = c("no decline",
                                paste0("decline ",
                                       sort(unique(decline[decline != 1]), decreasing = TRUE))))
  )

long <- res %>%
  select(Ne, Ne_num, Gen, decline, sel_s, af, n_parallel, f0, Ne_lab, Dec_lab,
         prob_any_false3, prob_any_false2, fdr3, fdr2) %>%
  pivot_longer(c(prob_any_false3, prob_any_false2, fdr3, fdr2),
               names_to = "key", values_to = "value") %>%
  mutate(
    overlap = ifelse(grepl("3$", key), "All three species", "At least two species"),
    overlap = factor(overlap, levels = names(overlap_cols)),
    metric  = ifelse(grepl("^prob", key), "P(at least one false positive)",
                     "False discovery rate"),
    metric  = factor(metric, levels = c("P(at least one false positive)",
                                        "False discovery rate"))
  )

s_vals <- sort(unique(res$sel_s))

fp_theme <- function(base_size, y_breaks = seq(0, 1, 0.25)) list(
  scale_colour_manual(values = overlap_cols, name = "Overlap"),
  scale_x_log10(breaks = s_vals, labels = as.character(s_vals)),
  scale_y_continuous(limits = c(0, 1), breaks = y_breaks),
  labs(x = "Selection coefficient (s, log scale)", y = NULL),
  theme_minimal(base_size = base_size),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
    strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
    strip.text.y     = element_text(angle = 0, hjust = 0),
    axis.text        = element_text(colour = "#6b6a66"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    legend.position  = "top"
  )
)

# ---------------- plots ----------------
for (n_par in sort(unique(res$n_parallel))) {
  for (g in sort(unique(res$Gen))) {

    od <- file.path(out_root, sprintf("par%03d", n_par), "false_positives", paste0("Gen_", g))
    dir.create(od, recursive = TRUE, showWarnings = FALSE)

    dg <- long %>% filter(n_parallel == n_par, Gen == g)

    # 1. Main: both metrics (rows) x f0 (columns), baseline Ne / decline
    d_main <- dg %>% filter(Ne == base_Ne, decline == base_decline)
    if (nrow(d_main) > 0) {
      p_main <- ggplot(d_main, aes(sel_s, value, colour = overlap, group = overlap)) +
        geom_hline(data = chance, aes(yintercept = value, colour = overlap),
                   linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
        geom_line(linewidth = 0.8, na.rm = TRUE) +
        geom_point(size = 1.8, na.rm = TRUE) +
        facet_grid(metric ~ f0, switch = "y") +
        fp_theme(size_main) +
        theme(strip.placement = "outside") +
        labs(title = sprintf("False positives in the cross-species overlap | %s generations", g),
             subtitle = sprintf(paste0("Ne = %s, no decline; %d parallel genes. ",
                                       "Dashed: chance with no selection. ",
                                       "FDR is blank where the overlap was always empty."),
                                base_Ne, n_par))
      ggsave(file.path(od, "fp_curves.main.png"), p_main,
             width = 14, height = 6.5, dpi = 300, bg = "white")
    }

    # 2. All conditions: P(>=1 FP), rows Ne > decline, columns f0
    d_all <- dg %>% filter(metric == "P(at least one false positive)")
    n_rows <- n_distinct(d_all$Ne_lab) * n_distinct(d_all$Dec_lab)
    p_all <- ggplot(d_all, aes(sel_s, value, colour = overlap, group = overlap)) +
      geom_hline(data = chance, aes(yintercept = value, colour = overlap),
                 linetype = "dashed", linewidth = 0.35, show.legend = FALSE) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.1) +
      facet_grid(Ne_lab + Dec_lab ~ f0) +
      fp_theme(size_all, y_breaks = c(0, 0.5, 1)) +
      labs(y = "P(at least one false positive)",
           title = sprintf("False positives, all conditions | %s generations", g),
           subtitle = "Dashed: chance with no selection")
    ggsave(file.path(od, "fp_curves.all_conditions.png"), p_all,
           width = 14, height = 1.7 * n_rows + 1.8, dpi = 300, bg = "white",
           limitsize = FALSE)

    message("Saved: ", od)
  }
}