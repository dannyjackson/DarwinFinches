#!/usr/bin/env Rscript
# ============================================================
# If a gene is detected in all three species, what is the probability that
# it is truly under parallel selection?
#
#   detected      = in the top-TOP_N FST list of all three species
#   true positive = selected in all three species (a "parallel" gene)
#
#   P(true | detected) = sum over permutations of parallel genes detected
#                        / sum over permutations of all genes detected
#
# Permutations where no gene was detected contribute nothing (0 / 0) and are
# ignored. Conditions where nothing was ever detected are left blank.
#
# Error bars: 95% Wilson interval treating each detected gene as one trial
# (approximate: genes in the same permutation are not fully independent).
# Point size: how many genes were detected in total, i.e. how much the
# estimate is based on.
#
# Input:  output/parallel_overlap.permutations.tsv  (parallel_outlier_overlap.R)
# Run from the project directory:
#   Rscript ${SCRIPTDIR}/plot_ppv_curves.R
# Writes:
#   output/parallel_overlap.ppv.tsv
#   output/power_curves_parallel/par<NN>/ppv/Gen_<g>/
#       ppv_curves.main.png             baseline Ne / no decline
#       ppv_curves.all_conditions.png   rows: Ne > decline, columns: f0
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

perm_file <- "output/parallel_overlap.permutations.tsv"
out_root  <- "output/power_curves_parallel"

base_Ne      <- "100k"
base_decline <- 1

size_main <- 13
size_all  <- 10

line_col <- "#2a78d6"

# ---------------- data ----------------
wilson <- function(x, n, z = 1.96) {
  p <- ifelse(n > 0, x / n, NA_real_)
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  list(lo = ifelse(n > 0, ctr - half, NA_real_),
       hi = ifelse(n > 0, ctr + half, NA_real_))
}

ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

perm <- read.delim(perm_file)
for (col in c("overlap3_n", "overlap3_parallel")) {
  if (!col %in% names(perm)) stop("Missing column ", col, " in ", perm_file)
}

ppv <- perm %>%
  group_by(Ne, Gen, decline, sel_s, af, n_parallel) %>%
  summarise(
    n_permutations          = n(),
    n_perm_with_detection   = sum(overlap3_n > 0),
    genes_detected          = sum(overlap3_n),
    true_detected           = sum(overlap3_parallel),
    false_detected          = sum(overlap3_n - overlap3_parallel),
    .groups = "drop"
  ) %>%
  mutate(
    ppv = ifelse(genes_detected > 0, true_detected / genes_detected, NA_real_),
    prob_false_given_detected = 1 - ppv
  )

ci <- wilson(ppv$true_detected, ppv$genes_detected)
ppv$ppv_lo <- ci$lo
ppv$ppv_hi <- ci$hi

write.table(ppv, "output/parallel_overlap.ppv.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: output/parallel_overlap.ppv.tsv")

cat(sprintf("Overall: %s of %s detected genes were true parallel genes (%.4f); %s false.\n",
            format(sum(ppv$true_detected), big.mark = ","),
            format(sum(ppv$genes_detected), big.mark = ","),
            sum(ppv$true_detected) / max(1, sum(ppv$genes_detected)),
            format(sum(ppv$false_detected), big.mark = ",")))

res <- ppv %>%
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

s_vals <- sort(unique(res$sel_s))

ppv_layers <- function(base_size, y_breaks = seq(0, 1, 0.25), pt_range = c(1.2, 4)) list(
  geom_errorbar(aes(ymin = ppv_lo, ymax = ppv_hi), width = 0.06,
                colour = line_col, alpha = 0.6, na.rm = TRUE),
  geom_line(colour = line_col, linewidth = 0.7, na.rm = TRUE),
  geom_point(aes(size = genes_detected), colour = line_col, na.rm = TRUE),
  scale_size_continuous(range = pt_range, name = "Genes detected\n(all permutations)"),
  scale_x_log10(breaks = s_vals, labels = as.character(s_vals)),
  scale_y_continuous(limits = c(0, 1), breaks = y_breaks),
  labs(x = "Selection coefficient (s, log scale)",
       y = "P(true parallel gene | detected in all three)"),
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

for (n_par in sort(unique(res$n_parallel))) {
  for (g in sort(unique(res$Gen))) {

    od <- file.path(out_root, sprintf("par%03d", n_par), "ppv", paste0("Gen_", g))
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
    dg <- res %>% filter(n_parallel == n_par, Gen == g)

    d_main <- dg %>% filter(Ne == base_Ne, decline == base_decline)
    if (nrow(d_main) > 0) {
      p_main <- ggplot(d_main, aes(sel_s, ppv, group = 1)) +
        ppv_layers(size_main) +
        facet_grid(. ~ f0) +
        labs(title = sprintf("If a gene is detected in all three species, is it real? | %s generations", g),
             subtitle = sprintf(paste0("Ne = %s, no decline; %d parallel genes. ",
                                       "Blank = no gene ever detected. Bars: 95%% interval."),
                                base_Ne, n_par))
      ggsave(file.path(od, "ppv_curves.main.png"), p_main,
             width = 14, height = 4.6, dpi = 300, bg = "white")
    }

    n_rows <- n_distinct(dg$Ne_lab) * n_distinct(dg$Dec_lab)
    p_all <- ggplot(dg, aes(sel_s, ppv, group = 1)) +
      ppv_layers(size_all, y_breaks = c(0, 0.5, 1), pt_range = c(0.6, 2.5)) +
      facet_grid(Ne_lab + Dec_lab ~ f0) +
      labs(title = sprintf("P(true parallel gene | detected), all conditions | %s generations", g),
           subtitle = "Blank = no gene ever detected. Bars: 95% interval.")
    ggsave(file.path(od, "ppv_curves.all_conditions.png"), p_all,
           width = 14, height = 1.7 * n_rows + 1.8, dpi = 300, bg = "white",
           limitsize = FALSE)

    message("Saved: ", od)
  }
}