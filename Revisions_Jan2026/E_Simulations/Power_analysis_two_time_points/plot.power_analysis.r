#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
root_dir <- "/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis"   # change if needed
out_dir  <- file.path(root_dir, "plots_fst_deltaTajD_neutral_vs_selected")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# optional filters (set to NULL for everything)
keep_ne      <- NULL  # e.g. c(10000, 100000)
keep_gen     <- NULL  # e.g. c(5, 15)
keep_decline <- NULL  # e.g. c(0.80, 1.00)
keep_f0      <- NULL  # only applies if a middle token exists & is numeric
keep_sel     <- NULL  # e.g. c(0.00, 0.01, 0.05, 0.10, 0.20)

# If you want to pin x-min (set to NA to not pin)
xmin_pin <- NA_real_

# Ellipses (set FALSE if you only want points)
draw_ellipses <- TRUE
ellipse_level <- 0.95
ellipse_min_n <- 10

# Point styling
pt_size  <- 0.55
pt_alpha <- 0.75

# ------------------------------------------------------------
# Discover files
# ------------------------------------------------------------
files <- Sys.glob(file.path(root_dir, "Ne_*", "Gen_*", "*.summary_stats.tsv"))
message("Found ", length(files), " summary_stats.tsv files")
if (length(files) == 0) stop("No *.summary_stats.tsv files found under: ", root_dir)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_ne <- function(ne_chr) {
  # accepts "Ne_10k", "Ne_10000", "Ne_100k"
  x <- str_remove(ne_chr, "^Ne_")
  if (str_detect(x, "k$")) as.numeric(str_remove(x, "k$")) * 1000 else suppressWarnings(as.numeric(x))
}

# Parse filename of form:
#   <decline>.<maybe_f0>.<sel>.summary_stats.tsv
# but also robust to:
#   <decline>.<sel>.summary_stats.tsv
# and extra dots; we take:
#   first numeric token as decline,
#   last numeric token as sel_s,
#   everything in between (if numeric) as f0 (joined with '.')
parse_from_filename <- function(fname) {
  base <- basename(fname)
  core <- str_remove(base, "\\.summary_stats\\.tsv$")

  toks <- str_split(core, "\\.", simplify = TRUE)
  toks <- toks[1, toks[1, ] != ""]
  if (length(toks) < 2) return(list(decline = NA_real_, f0 = NA_real_, sel_s = NA_real_))

  # decline: first two tokens joined if that yields something like 0.80 / 1.00
  # fallback: first token numeric
  decline_try <- suppressWarnings(as.numeric(paste0(toks[1], ".", toks[2])))
  if (is.finite(decline_try)) {
    decline <- decline_try
    # remove the two tokens used
    rem <- toks[-c(1, 2)]
  } else {
    decline <- suppressWarnings(as.numeric(toks[1]))
    rem <- toks[-1]
  }

  # sel: last token (or last two joined if needed)
  sel <- suppressWarnings(as.numeric(rem[length(rem)]))
  rem_mid <- rem[-length(rem)]

  # middle: if anything remains, join and parse as numeric (often 0 / 0.0 / 0.05 etc.)
  f0 <- NA_real_
  if (length(rem_mid) > 0) {
    mid_str <- paste(rem_mid, collapse = ".")
    mid_num <- suppressWarnings(as.numeric(mid_str))
    if (is.finite(mid_num)) f0 <- mid_num
  }

  list(decline = decline, f0 = f0, sel_s = sel)
}

# ------------------------------------------------------------
# Parse metadata from path
# ------------------------------------------------------------
parse_from_path <- function(path) {
  rel <- gsub(paste0("^", root_dir, "/"), "", path)
  parts <- str_split(rel, "/", simplify = TRUE)

  # expected:
  # [1] Ne_100k
  # [2] Gen_15
  # [3] 0.80.0.00.summary_stats.tsv
  if (ncol(parts) < 3) return(tibble())

  ne_chr  <- parts[1]
  gen_chr <- parts[2]
  fn      <- parts[3]

  ne_val  <- parse_ne(ne_chr)
  gen_val <- suppressWarnings(as.numeric(str_remove(gen_chr, "^Gen_")))

  fp <- parse_from_filename(fn)

  tibble(
    file     = path,
    Ne       = ne_val,
    Gen_time = gen_val,
    decline  = fp$decline,
    f0       = fp$f0,
    sel_s    = fp$sel_s
  )
}

meta <- map_dfr(files, parse_from_path) %>%
  filter(!is.na(Ne), !is.na(Gen_time), !is.na(decline), !is.na(sel_s))

if (nrow(meta) == 0) {
  stop("Parsed 0 metadata rows. Check that files match Ne_*/Gen_*/*.summary_stats.tsv and filenames encode decline/sel.")
}

message("Parsed metadata for ", nrow(meta), " files")
message("Unique parameter values:")
message("  Ne:        ", paste(sort(unique(meta$Ne)), collapse = ", "))
message("  Gen_time:  ", paste(sort(unique(meta$Gen_time)), collapse = ", "))
message("  decline:   ", paste(sort(unique(meta$decline)), collapse = ", "))
message("  f0:        ", paste(sort(unique(meta$f0[is.finite(meta$f0)])), collapse = ", "))
message("  sel_s:     ", paste(sort(unique(meta$sel_s)), collapse = ", "))

# Optional filters
if (!is.null(keep_ne))      meta <- meta %>% filter(Ne %in% keep_ne)
if (!is.null(keep_gen))     meta <- meta %>% filter(Gen_time %in% keep_gen)
if (!is.null(keep_decline)) meta <- meta %>% filter(decline %in% keep_decline)
if (!is.null(keep_sel))     meta <- meta %>% filter(sel_s %in% keep_sel)
if (!is.null(keep_f0))      meta <- meta %>% filter(is.na(f0) | f0 %in% keep_f0)

if (nrow(meta) == 0) stop("After filtering, no files remain.")

# ------------------------------------------------------------
# Robust reader (drops malformed lines by field count)
# ------------------------------------------------------------
forced_types <- cols(
  rep           = col_character(),
  sel_s         = col_double(),
  decline_rate  = col_double(),
  offset        = col_double(),
  model         = col_character(),
  seed          = col_character(),
  num_sites     = col_double(),
  num_mutations = col_double(),
  pi_t1         = col_double(),
  pi_t2         = col_double(),
  dxy           = col_double(),
  fst_hudson    = col_double(),
  tajd_t1       = col_double(),
  tajd_t2       = col_double(),
  delta_tajd    = col_double(),
  .default      = col_guess()
)

read_one <- function(file, Ne, Gen_time, decline, f0, sel_s) {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) < 2) return(tibble())

  expected_n <- length(strsplit(lines[1], "\t", fixed = TRUE)[[1]])
  good <- vapply(strsplit(lines, "\t", fixed = TRUE), length, integer(1)) == expected_n

  df <- readr::read_tsv(I(lines[good]), col_types = forced_types, show_col_types = FALSE)

  need <- c("fst_hudson", "delta_tajd")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) stop("Missing required columns in ", file, ": ", paste(missing, collapse = ", "))

  df %>%
    mutate(
      file     = file,
      Ne       = Ne,
      Gen_time = Gen_time,
      decline  = decline,
      f0       = f0,
      sel_s    = sel_s
    )
}

df_all <- pmap_dfr(meta, read_one)
message("Finished reading data. Total rows loaded: ", nrow(df_all))
if (nrow(df_all) == 0) stop("No rows read from files.")

# ------------------------------------------------------------
# Make neutral-vs-selected label
# ------------------------------------------------------------
df_all <- df_all %>%
  mutate(
    sel_num = as.numeric(sel_s),
    sel_f   = factor(sprintf("%.2f", sel_num),
                     levels = sprintf("%.2f", sort(unique(sel_num))))
  )

# ------------------------------------------------------------
# Color map (must happen AFTER sel_f exists)
# ------------------------------------------------------------
sel_levels <- levels(df_all$sel_f)

sel_cols <- c(
  "0.00" = "grey60",
  "0.01" = "#F4A6A6",
  "0.05" = "#F26B6B",
  "0.10" = "#D94A4A",
  "0.20" = "#8B0000"
)

# fallback for unexpected s values:
missing_levels <- setdiff(sel_levels, names(sel_cols))
if (length(missing_levels) > 0) sel_cols[missing_levels] <- "#B22222"

# ------------------------------------------------------------
# Plot per (Ne, Gen_time, decline, f0)
# ------------------------------------------------------------
# If f0 is entirely NA, it won't create separate groups by f0.
grp_vars <- c("Ne", "Gen_time", "decline", "f0")

df_split <- df_all %>%
  group_by(across(all_of(grp_vars))) %>%
  group_split()

message("Preparing plots for ", length(df_split), " parameter combinations")

make_plot <- function(d) {
  # ensure both neutral and selected are present; still plot if only one exists
  ellipse_df <- d %>%
    group_by(sel_f) %>%
    filter(n() >= ellipse_min_n) %>%
    ungroup()

  p <- ggplot(d, aes(x = fst_hudson, y = delta_tajd)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70")

  if (draw_ellipses && nrow(ellipse_df) > 0) {
  p <- p +
    stat_ellipse(
      data = ellipse_df,
      aes(color = sel_f, group = sel_f),
      type  = "norm",
      level = ellipse_level,
      linewidth = 0.45,
      alpha = 0.9
    )
}

p <- p +
  geom_point(
    aes(color = sel_f),
    size = pt_size,
    alpha = pt_alpha
  ) +
  scale_color_manual(values = sel_cols, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(color = "Selection (s)") +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom"
    ) +
    labs(
      x = "FST (Hudson)",
      y = expression(Delta * " Tajima's D"),
      color = NULL
    )

  if (!is.na(xmin_pin)) p <- p + coord_cartesian(xlim = c(xmin_pin, NA))
  p
}

for (d in df_split) {
  if (nrow(d) == 0) next

  Ne_val      <- unique(d$Ne)
  Gen_val     <- unique(d$Gen_time)
  decline_val <- unique(d$decline)
  f0_val      <- unique(d$f0)

  # Handle potentially NA f0
  f0_tag <- if (all(is.na(f0_val))) {
    "f0NA"
  } else {
    paste0("f0", formatC(f0_val[1], format = "f", digits = 2))
  }

  p <- make_plot(d) +
    ggtitle(paste0(
      "FST vs ΔTajima’s D | Ne=", Ne_val,
      " | Gen_time=", Gen_val,
      " | decline=", sprintf("%.2f", decline_val[1]),
      if (!all(is.na(f0_val))) paste0(" | f0=", sprintf("%.2f", f0_val[1])) else ""
    ))

  out_file <- file.path(
    out_dir,
    paste0(
      "fst_vs_deltaTajd_neutral_vs_selected_",
      "Ne", Ne_val,
      "_Gen", Gen_val,
      "_decl", sprintf("%.2f", decline_val[1]),
      "_", f0_tag,
      ".png"
    )
  )

  ggsave(filename = out_file, plot = p, width = 7.5, height = 6.0, dpi = 300)
  message("Saved: ", out_file)
}

message("\nDone. Wrote plots to: ", out_dir)




##################################################################
################# Make table of power analyses ###################
##################################################################
# ------------------------------------------------------------
# Outlier enrichment / "power": among top X% score, what fraction are selected?
# Run three identical analyses:
#   1) composite: z(FST) + -z(ΔTajD)
#   2) fst_only : z(FST)
#   3) taj_only : -z(ΔTajD)
# z computed within (neutral + one s) sample
# ------------------------------------------------------------

fst_col <- "fst_hudson"
taj_col <- "delta_tajd"     # or "tajd_t2"
top_q   <- 0.99             # NOTE: 0.99 = top 1%. For top 0.1%, use 0.999

# ensure numeric selection coefficient
df_all <- df_all %>%
  mutate(sel_num = as.numeric(as.character(sel_s)))

group_vars <- c("Ne", "Gen_time", "decline", "f0")

compute_enrichment <- function(df, score_type = c("composite", "fst_only", "taj_only"),
                               fst_col = "fst_hudson", taj_col = "delta_tajd", top_q = 0.99,
                               group_vars = c("Ne","Gen_time","decline","f0")) {
  score_type <- match.arg(score_type)

  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~{
      d <- .x

      # must have neutral
      if (!any(d$sel_num == 0)) return(tibble())

      s_vals <- sort(unique(d$sel_num[d$sel_num > 0]))
      if (length(s_vals) == 0) return(tibble())

      bind_rows(lapply(s_vals, function(sv) {
        dd <- d %>% filter(sel_num == 0 | sel_num == sv)

        z_fst <- as.numeric(scale(dd[[fst_col]]))
        z_taj <- as.numeric(scale(dd[[taj_col]]))

        score <- switch(
          score_type,
          composite = z_fst + (-1 * z_taj),
          fst_only  = z_fst,
          taj_only  = (-1 * z_taj)
        )

        cutoff <- as.numeric(quantile(score, probs = top_q, na.rm = TRUE))
        is_out <- score >= cutoff

        n_out     <- sum(is_out, na.rm = TRUE)
        n_sel_out <- sum(is_out & dd$sel_num == sv, na.rm = TRUE)

        tibble(
          score_type = score_type,
          sel_s = sv,
          n_total = nrow(dd),
          n_outliers = n_out,
          n_selected_outliers = n_sel_out,
          freq_selected = ifelse(n_out > 0, n_sel_out / n_out, NA_real_)
        )
      }))
    }) %>%
    ungroup() %>%
    arrange(Ne, Gen_time, decline, f0, sel_s)
}

# ---- Run all three ----
enrichment_composite <- compute_enrichment(df_all, "composite", fst_col, taj_col, top_q, group_vars)
enrichment_fst_only  <- compute_enrichment(df_all, "fst_only",  fst_col, taj_col, top_q, group_vars)
enrichment_taj_only  <- compute_enrichment(df_all, "taj_only",  fst_col, taj_col, top_q, group_vars)

# ---- Write next to plots ----
out_tbl_comp <- file.path(out_dir, "freq_selected_among_top_composite.tsv")
out_tbl_fst  <- file.path(out_dir, "freq_selected_among_top_fst_only.tsv")
out_tbl_taj  <- file.path(out_dir, "freq_selected_among_top_deltaTajD_only.tsv")

readr::write_tsv(enrichment_composite, out_tbl_comp)
readr::write_tsv(enrichment_fst_only,  out_tbl_fst)
readr::write_tsv(enrichment_taj_only,  out_tbl_taj)

message("Wrote enrichment tables:")
message("  composite: ", out_tbl_comp)
message("  fst_only : ", out_tbl_fst)
message("  taj_only : ", out_tbl_taj)

print(head(enrichment_composite, 20))
print(head(enrichment_fst_only, 20))
print(head(enrichment_taj_only, 20))

# ------------------------------------------------------------
# Comparison of "power" across score types:
# power := freq_selected = P(selected | window is in top tail)
# Outputs:
#   - power_comparison_wide.tsv  (per group/sel_s)
#   - power_summary_by_sel.tsv   (overall mean power by sel_s)
# ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(readr)

# (Assumes you already created:
#  enrichment_composite, enrichment_fst_only, enrichment_taj_only)

# Stack long
enrichment_all <- bind_rows(
  enrichment_composite,
  enrichment_fst_only,
  enrichment_taj_only
)

id_vars <- c(group_vars, "sel_s")  # group_vars defined earlier

# Wide comparison table
power_cmp_wide <- enrichment_all %>%
  select(all_of(id_vars), score_type, n_outliers, n_selected_outliers, freq_selected) %>%
  pivot_wider(
    names_from = score_type,
    values_from = c(n_outliers, n_selected_outliers, freq_selected),
    names_glue = "{.value}_{score_type}"
  ) %>%
  # nicer column names (optional)
  rename(
    power_composite = freq_selected_composite,
    power_fst       = freq_selected_fst_only,
    power_taj       = freq_selected_taj_only,
    n_out_composite = n_outliers_composite,
    n_out_fst       = n_outliers_fst_only,
    n_out_taj       = n_outliers_taj_only,
    n_sel_out_composite = n_selected_outliers_composite,
    n_sel_out_fst       = n_selected_outliers_fst_only,
    n_sel_out_taj       = n_selected_outliers_taj_only
  ) %>%
  arrange(Ne, Gen_time, decline, f0, sel_s)

out_cmp <- file.path(out_dir, "power_comparison_wide.tsv")
write_tsv(power_cmp_wide, out_cmp)
message("Wrote power comparison table: ", out_cmp)
print(head(power_cmp_wide, 20))

# Optional: overall summary by sel_s and score_type
power_summary_by_sel <- enrichment_all %>%
  group_by(sel_s, score_type) %>%
  summarise(
    n_groups = n(),
    mean_power = mean(freq_selected, na.rm = TRUE),
    median_power = median(freq_selected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(sel_s, score_type)

out_sum <- file.path(out_dir, "power_summary_by_sel.tsv")
write_tsv(power_summary_by_sel, out_sum)
message("Wrote power summary table: ", out_sum)
print(power_summary_by_sel)