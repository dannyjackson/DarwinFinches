#!/usr/bin/env Rscript
# ============================================================
# Summarise and plot gsea_power_empirical.R output
#
# Run from the gsea_power_analysis directory:
#   Rscript ${SCRIPTDIR}/summarize_gsea_power_empirical.R
#
# Reads:  gsea_empirical_power/Ne_*/Gen_*/*.iterations.tsv
# Writes: gsea_empirical_power/summary.tsv
#         gsea_empirical_power/plots/<metric>.Gen_<g>.png
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

in_dir   <- "gsea_empirical_power"
plot_dir <- file.path(in_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "\\.iterations\\.tsv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No *.iterations.tsv files under ", in_dir)

d <- rbindlist(lapply(files, fread), fill = TRUE)

# Keep only complete iterations (cross-species row present)
d <- d[d[species == "cross_species", .(Ne, Gen, decline, sel_s, f0, term, iteration)],
       on = .(Ne, Gen, decline, sel_s, f0, term, iteration)]

keys <- c("term", "Ne", "Gen", "decline", "sel_s", "f0")

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

cross <- d[species == "cross_species", .(
  n_iter              = .N,
  median_pct_cross    = median(pct_NES, na.rm = TRUE),
  power_cross_top1    = mean(pct_NES <= 1, na.rm = TRUE),
  power_cross_top5    = mean(pct_NES <= 5, na.rm = TRUE),
  power_cross_top10   = mean(pct_NES <= 10, na.rm = TRUE),
  power_sig_all3      = mean(n_species_sig == 3),
  power_sig_any       = mean(n_species_sig >= 1),
  prop_positive_all3  = mean(n_species_positive == 3)
), by = keys]

per_sp <- dcast(
  d[species != "cross_species", .(
    power_sig  = mean(padj < 0.05, na.rm = TRUE),
    median_pct = median(pct_NES, na.rm = TRUE)
  ), by = c(keys, "species")],
  as.formula(paste(paste(keys, collapse = " + "), "~ species")),
  value.var = c("power_sig", "median_pct")
)

summary_tab <- merge(cross, per_sp, by = keys)
setorderv(summary_tab, keys)
fwrite(summary_tab, file.path(in_dir, "summary.tsv"), sep = "\t")
message("Wrote: ", file.path(in_dir, "summary.tsv"))
print(summary_tab[, .(Ne, Gen, decline, sel_s, f0, n_iter,
                      median_pct_cross, power_cross_top5, power_sig_all3)])

# Null control (s = 0) values, used as the reference line when present
null_ref <- summary_tab[sel_s == 0]
sel_tab  <- summary_tab[sel_s > 0]
if (nrow(sel_tab) == 0) stop("No selected (s > 0) conditions to plot yet.")

# ------------------------------------------------------------
# Power curves: x = s (log), facets = f0 (rows: Ne x decline)
# Same look as plot_power_curves.R
# ------------------------------------------------------------

s_ramp <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")

metrics <- data.table(
  column = c("power_cross_top5", "power_cross_top1", "power_sig_all3", "median_pct_cross"),
  label  = c("P(term in top 5% by mean NES across species)",
             "P(term in top 1% by mean NES across species)",
             "P(padj < 0.05 in all three species)",
             "Median percentile of term (top x%)"),
  expect = c(0.05, 0.01, NA, 50)     # expectation with no selection
)

ne_num <- function(x) ifelse(grepl("k$", x), as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

sel_tab[, `:=`(
  s       = factor(sel_s, levels = sort(unique(sel_s))),
  f0_lab  = factor(paste0("f0 = ", f0), levels = paste0("f0 = ", sort(unique(f0)))),
  Ne_lab  = paste0("Ne ", Ne),
  Dec_lab = ifelse(decline == 1, "no decline", paste0("decline ", decline))
)]
sel_tab[, Ne_lab := factor(Ne_lab, levels = unique(Ne_lab[order(ne_num(Ne))]))]
sel_tab[, Dec_lab := factor(Dec_lab, levels = c("no decline", sort(setdiff(unique(Dec_lab), "no decline"), decreasing = TRUE)))]

s_cols <- setNames(rep_len(s_ramp, nlevels(sel_tab$s)), levels(sel_tab$s))
s_vals <- sort(unique(sel_tab$sel_s))

for (m in seq_len(nrow(metrics))) {
  col <- metrics$column[m]
  for (g in sort(unique(sel_tab$Gen))) {
    dg <- sel_tab[Gen == g]
    n_rows <- uniqueN(dg$Ne_lab) * uniqueN(dg$Dec_lab)

    # Reference line: observed null control if it exists for this Gen, else expectation
    ref <- null_ref[Gen == g, mean(get(col), na.rm = TRUE)]
    ref_txt <- "expectation with no selection"
    if (is.na(ref) || !is.finite(ref)) ref <- metrics$expect[m] else ref_txt <- "observed in s = 0 control"

    p <- ggplot(dg, aes(sel_s, .data[[col]], group = 1)) +
      { if (!is.na(ref)) geom_hline(yintercept = ref, linetype = "dashed",
                                    colour = "#8a8984", linewidth = 0.35) } +
      geom_line(linewidth = 0.7, colour = "#3d3d3a") +
      geom_point(aes(colour = s), size = 1.8) +
      scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE) +
      scale_x_log10(breaks = s_vals, labels = as.character(s_vals)) +
      { if (col == "median_pct_cross")
          scale_y_reverse(limits = c(100, 0), breaks = c(0, 25, 50, 75, 100))
        else scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) } +
      facet_grid(Ne_lab + Dec_lab ~ f0_lab) +
      labs(x = "Selection coefficient (s, log scale)", y = metrics$label[m],
           title = sprintf("Power of empirical GSEA to detect %s | %s generations",
                           unique(dg$term)[1], g),
           subtitle = if (!is.na(ref)) sprintf("Dashed line: %s (%g)", ref_txt, signif(ref, 2)) else NULL) +
      theme_minimal(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
        strip.text.y     = element_text(angle = 0, hjust = 0),
        axis.text        = element_text(colour = "#6b6a66"),
        axis.text.x      = element_text(angle = 45, hjust = 1),
        plot.title       = element_text(face = "bold")
      )

    out <- file.path(plot_dir, sprintf("%s.Gen_%s.png", col, g))
    ggsave(out, p, width = 14, height = 2.2 * n_rows + 1.6, dpi = 300, bg = "white",
           limitsize = FALSE)
    message("Saved: ", out)
  }
}
