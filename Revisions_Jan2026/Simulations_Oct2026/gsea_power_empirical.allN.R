#!/usr/bin/env Rscript
# ============================================================
# Power of the empirical cross-species GSEA (blood coagulation) as a function
# of N, the number of genes in the term under selection.
#
# Reads:  gsea_empirical_power_nsel/Ne_*/Gen_*/*.nsel.iterations.tsv
#         (from gsea_power_empirical.R with the n_sel option)
# Writes: gsea_empirical_power_nsel/summary_nsel.tsv
#         gsea_empirical_power_nsel/plots/<metric>/Gen_<g>/
#             power_vs_nsel.main.png            baseline Ne / no decline
#             power_vs_nsel.all_conditions.png  rows: Ne > decline
#
# Each plot: x = number of selected genes in the term (N), y = power,
# one line per s, panels = f0. Dashed line = power with no selected genes
# (observed at N = 0, which is the null for every condition).
#
# Run from the project directory:
#   Rscript ${SCRIPTDIR}/plot_gsea_parallel_nsel.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_dir   <- "gsea_empirical_power_nsel"
plot_dir <- file.path(in_dir, "plots")

base_Ne      <- "100k"
base_decline <- 1

size_main <- 13
size_all  <- 10

# light = weak selection, dark = strong selection (up to 9 s values)
s_ramp <- c("#b7d3f6", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
            "#2a78d6", "#1c5cab", "#104281", "#0d366b")

# ---------------- read and summarise ----------------
files <- list.files(in_dir, pattern = "\\.nsel\\.iterations\\.tsv$",
                    recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No *.nsel.iterations.tsv files under ", in_dir)

d <- rbindlist(lapply(files, fread), fill = TRUE)
keys <- c("term", "Ne", "Gen", "decline", "sel_s", "f0", "n_selected")

# keep only complete (iteration, N) pairs: cross-species row present
d <- d[d[species == "cross_species", .(Ne, Gen, decline, sel_s, f0, term, n_selected, iteration)],
       on = .(Ne, Gen, decline, sel_s, f0, term, n_selected, iteration)]

cross <- d[species == "cross_species", .(
  n_iter            = .N,
  median_pct_cross  = median(pct_NES, na.rm = TRUE),
  power_cross_top1  = mean(pct_NES <= 1, na.rm = TRUE),
  power_cross_top5  = mean(pct_NES <= 5, na.rm = TRUE),
  power_cross_top10 = mean(pct_NES <= 10, na.rm = TRUE),
  power_sig_all3    = mean(n_species_sig == 3),
  power_sig_any     = mean(n_species_sig >= 1)
), by = keys]

per_sp <- dcast(
  d[species != "cross_species", .(power_sig = mean(padj < 0.05, na.rm = TRUE),
                                  median_pct = median(pct_NES, na.rm = TRUE)),
    by = c(keys, "species")],
  as.formula(paste(paste(keys, collapse = " + "), "~ species")),
  value.var = c("power_sig", "median_pct")
)

summ <- merge(cross, per_sp, by = keys)
setorderv(summ, keys)
fwrite(summ, file.path(in_dir, "summary_nsel.tsv"), sep = "\t")
message("Wrote: ", file.path(in_dir, "summary_nsel.tsv"))

cat("\nN values found:", paste(sort(unique(summ$n_selected)), collapse = ", "), "\n")
cat("Iterations per point: min", min(summ$n_iter), "max", max(summ$n_iter), "\n")

# ---------------- null reference (N = 0) ----------------
metrics <- data.table(
  column = c("power_cross_top5", "power_cross_top1", "power_sig_all3"),
  tag    = c("cross_top5",       "cross_top1",       "sig_all3"),
  label  = c("P(term in top 5% by mean NES across species)",
             "P(term in top 1% by mean NES across species)",
             "P(padj < 0.05 in all three species)"),
  expect = c(0.05, 0.01, 0)
)

null_tab <- summ[n_selected == 0]
cat("\nObserved power with no selected genes (N = 0), by generation time:\n")
print(null_tab[, lapply(.SD, function(x) round(weighted.mean(x, n_iter), 3)),
               by = Gen, .SDcols = metrics$column])

# ---------------- plotting data ----------------
ne_num <- function(x) ifelse(grepl("k$", x), as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

# Selected conditions only (s > 0); N = 0 points kept on each curve as its start
sel <- summ[sel_s > 0]
if (nrow(sel) == 0) stop("No selected (s > 0) conditions yet.")
sel[, `:=`(
  s       = factor(sel_s, levels = sort(unique(sel_s))),
  f0_lab  = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", sort(unique(f0)))),
  Ne_lab  = paste0("Ne ", Ne),
  Dec_lab = ifelse(decline == 1, "no decline", paste0("decline ", decline))
)]
sel[, Ne_lab := factor(Ne_lab, levels = unique(Ne_lab[order(ne_num(Ne))]))]
sel[, Dec_lab := factor(Dec_lab, levels = c("no decline",
                                            sort(setdiff(unique(Dec_lab), "no decline"),
                                                 decreasing = TRUE)))]

s_cols <- setNames(rep_len(s_ramp[seq(length(s_ramp) - nlevels(sel$s) + 1, length(s_ramp))],
                           nlevels(sel$s)), levels(sel$s))
n_max  <- max(sel$n_selected)
term   <- unique(sel$term)[1]

theme_nsel <- function(base_size) theme_minimal(base_size = base_size) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
        strip.text.y     = element_text(angle = 0, hjust = 0),
        axis.text        = element_text(colour = "#6b6a66"),
        plot.title       = element_text(face = "bold"),
        legend.position  = "right")

for (m in seq_len(nrow(metrics))) {
  col <- metrics$column[m]

  for (g in sort(unique(sel$Gen))) {
    od <- file.path(plot_dir, metrics$tag[m], paste0("Gen_", g))
    dir.create(od, recursive = TRUE, showWarnings = FALSE)

    ref <- null_tab[Gen == g, weighted.mean(get(col), n_iter)]
    ref_txt <- "observed with no selected genes (N = 0)"
    if (length(ref) == 0 || is.na(ref)) { ref <- metrics$expect[m]; ref_txt <- "expectation with no selection" }

    dg <- sel[Gen == g]

    base_layers <- list(
      geom_hline(yintercept = ref, linetype = "dashed", colour = "#8a8984", linewidth = 0.4),
      geom_line(linewidth = 0.8),
      geom_point(size = 1.6),
      scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE),
      scale_x_continuous(breaks = sort(unique(sel$n_selected)), limits = c(0, n_max)),
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)),
      labs(x = sprintf("Genes in %s under selection (N, of %d)", term, n_max),
           y = metrics$label[m])
    )

    d_main <- dg[Ne == base_Ne & abs(decline - base_decline) < 1e-9]
    if (nrow(d_main) > 0) {
      p <- ggplot(d_main, aes(n_selected, .data[[col]], colour = s, group = s)) +
        base_layers + facet_grid(. ~ f0_lab) + theme_nsel(size_main) +
        labs(title = sprintf("Cross-species GSEA power vs. number of selected genes | %s generations", g),
             subtitle = sprintf("Ne = %s, no decline. Dashed line: %s (%.3f).",
                                base_Ne, ref_txt, ref))
      ggsave(file.path(od, "power_vs_nsel.main.png"), p,
             width = 15, height = 4.6, dpi = 300, bg = "white")
    }

    n_rows <- uniqueN(dg$Ne_lab) * uniqueN(dg$Dec_lab)
    p_all <- ggplot(dg, aes(n_selected, .data[[col]], colour = s, group = s)) +
      base_layers + facet_grid(Ne_lab + Dec_lab ~ f0_lab) + theme_nsel(size_all) +
      labs(title = sprintf("Cross-species GSEA power vs. N, all conditions | %s generations", g),
           subtitle = sprintf("Dashed line: %s (%.3f)", ref_txt, ref))
    ggsave(file.path(od, "power_vs_nsel.all_conditions.png"), p_all,
           width = 15, height = 1.9 * n_rows + 1.8, dpi = 300, bg = "white",
           limitsize = FALSE)

    message("Saved: ", od)
  }
}