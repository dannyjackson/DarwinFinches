#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
})

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------
root_dir <- "output"
out_dir  <- file.path(root_dir, "plots_fst_deltaTajD")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# optional filters (set to NULL for everything)
keep_f0      <- NULL  # e.g. c(0.05, 0.25)
keep_decline <- NULL  # e.g. c(0.00, 0.10, 0.20)
keep_sel     <- NULL  # e.g. c(0.00, 0.01, 0.05, 0.10, 0.20)
keep_ne      <- NULL  # e.g. c(10000, 50000)
keep_gen     <- NULL  # e.g. c(5, 10, 15)

# If you want to pin the xmin, change this (or set to NA to not pin)
xmin_pin <- -0.05

# ------------------------------------------------------------
# Find all summary_stats.tsv files
# Expected layout:
# output/Ne_10k/Gen_5/AF_0_05/0.00/selection_0.10/summary_stats.tsv
# ------------------------------------------------------------
files <- Sys.glob(file.path(
  root_dir,
  "Ne_*", "Gen_*", "AF_*", "*", "selection_*", "summary_stats.tsv"
))


message("Found ", length(files), " summary_stats.tsv files")

if (length(files) == 0) stop("No summary_stats.tsv files found under: ", root_dir)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
parse_ne <- function(ne_chr) {
  # accepts "Ne_10k", "Ne_10000"
  x <- str_remove(ne_chr, "^Ne_")
  if (str_detect(x, "k$")) {
    as.numeric(str_remove(x, "k$")) * 1000
  } else {
    suppressWarnings(as.numeric(x))
  }
}

parse_f0 <- function(af_chr) {
  # "AF_0_05" -> 0.05, "AF_0_95" -> 0.95
  x <- str_remove(af_chr, "^AF_")
  x <- str_replace_all(x, "_", ".")
  suppressWarnings(as.numeric(x))
}

# ------------------------------------------------------------
# Parse metadata from path
# ------------------------------------------------------------
parse_from_path <- function(path) {
  rel <- gsub(paste0("^", root_dir, "/"), "", path)
  parts <- str_split(rel, "/", simplify = TRUE)

  # parts should look like:
  # [1] Ne_10k
  # [2] Gen_5
  # [3] AF_0_05
  # [4] 0.00              <-- decline dir
  # [5] selection_0.10
  # [6] summary_stats.tsv

  if (ncol(parts) < 6) return(tibble())

  ne_chr      <- parts[1]
  gen_chr     <- parts[2]
  af_chr      <- parts[3]
  decline_chr <- parts[4]

  sel_chr <- str_match(rel, "/selection_([0-9]+\\.[0-9]+|[0-9]+)")[, 2]

  tibble(
    file        = path,
    Ne          = parse_ne(ne_chr),
    Gen_time    = suppressWarnings(as.numeric(str_remove(gen_chr, "^Gen_"))),
    f0          = parse_f0(af_chr),
    decline     = suppressWarnings(as.numeric(decline_chr)),
    sel_s       = suppressWarnings(as.numeric(sel_chr))
  )
}

meta <- map_dfr(files, parse_from_path) %>%
  filter(!is.na(Ne), !is.na(Gen_time), !is.na(f0), !is.na(decline), !is.na(sel_s))

message("Parsed metadata for ", nrow(meta), " files")

message("Unique parameter values:")
message("  Ne:        ", paste(sort(unique(meta$Ne)), collapse = ", "))
message("  Gen_time:  ", paste(sort(unique(meta$Gen_time)), collapse = ", "))
message("  f0:        ", paste(sort(unique(meta$f0)), collapse = ", "))
message("  decline:   ", paste(sort(unique(meta$decline)), collapse = ", "))
message("  sel_s:     ", paste(sort(unique(meta$sel_s)), collapse = ", "))

if (nrow(meta) == 0) {
  stop(
    "Found summary_stats.tsv files, but none matched expected path structure:\n",
    "  output/Ne_*/Gen_*/AF_*/<decline>/selection_*/summary_stats.tsv"
  )
}

# Optional filters
if (!is.null(keep_f0))      meta <- meta %>% filter(f0 %in% keep_f0)
if (!is.null(keep_decline)) meta <- meta %>% filter(decline %in% keep_decline)
if (!is.null(keep_sel))     meta <- meta %>% filter(sel_s %in% keep_sel)
if (!is.null(keep_ne))      meta <- meta %>% filter(Ne %in% keep_ne)
if (!is.null(keep_gen))     meta <- meta %>% filter(Gen_time %in% keep_gen)

if (nrow(meta) == 0) stop("After filtering, no files remain to plot.")

# ------------------------------------------------------------
# Reader with robust row filtering based on header field count
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

read_one <- function(file, Ne, Gen_time, f0, decline, sel_s) {
  lines <- readLines(file, warn = FALSE)
  if (length(lines) < 2) return(tibble())

  # infer expected number of tab-fields from header
  expected_n <- length(strsplit(lines[1], "\t", fixed = TRUE)[[1]])
  good <- vapply(strsplit(lines, "\t", fixed = TRUE), length, integer(1)) == expected_n

  df <- readr::read_tsv(I(lines[good]), col_types = forced_types, show_col_types = FALSE)

  need <- c("fst_hudson", "delta_tajd", "sel_s")
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) stop("Missing required columns in ", file, ": ", paste(missing, collapse = ", "))

  df %>%
    mutate(
      file     = file,
      Ne       = Ne,
      Gen_time = Gen_time,
      f0       = f0,
      decline  = decline,
      sel_s    = sel_s
    )
}

df_all <- pmap_dfr(meta, read_one)

message("Finished reading data")
message("Total rows loaded: ", nrow(df_all))

if (nrow(df_all) == 0) stop("No rows read from files (all empty or filtered out).")

# ------------------------------------------------------------
# Factor ordering for facet grid (rows=f0, cols=decline)
# ------------------------------------------------------------
f0_levels  <- sort(unique(df_all$f0), decreasing = FALSE)
dec_levels <- sort(unique(df_all$decline), decreasing = FALSE)

# consistent selection ordering + fixed labels
sel_levels <- sort(unique(df_all$sel_s), decreasing = FALSE)

df_all <- df_all %>%
  mutate(
    f0      = factor(f0, levels = f0_levels, labels = sprintf("%.2f", f0_levels)),
    decline = factor(decline, levels = dec_levels, labels = sprintf("%.2f", dec_levels)),
    sel_s   = factor(sel_s, levels = sel_levels, labels = sprintf("%.2f", sel_levels))
  )

# ------------------------------------------------------------
# Color map: enforce your requested colors, fill in the others
# ------------------------------------------------------------
sel_cols <- c(
  "0.00" = "grey60",
  "0.01" = "#F4A6A6",  # light red
  "0.05" = "#F26B6B",
  "0.10" = "#D94A4A",
  "0.20" = "#8B0000"   # dark red
)

auc_rank <- function(scores, labels01) {
  # labels01: 0/1
  r <- rank(scores, ties.method = "average")
  n1 <- sum(labels01 == 1)
  n0 <- sum(labels01 == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[labels01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# If some sel values exist beyond these, give them a fallback:
missing_levels <- setdiff(levels(df_all$sel_s), names(sel_cols))
if (length(missing_levels) > 0) {
  # fallback to a mid red for any unexpected levels
  sel_cols[missing_levels] <- "#B22222"
}

# ------------------------------------------------------------
# Plot-per-(Ne, Gen_time)
# ------------------------------------------------------------

df_split <- split(df_all, list(df_all$Ne, df_all$Gen_time), drop = TRUE)

message("Preparing plots for ", length(df_split), " (Ne × Gen_time) combinations")

make_plot <- function(d) {
  # numeric sel for logic
  d <- d %>%
    mutate(sel_num = as.numeric(as.character(sel_s)))

    # ---- separation metric per facet: AUC of Mahalanobis-to-neutral
    # ---- Per-s separation metric per facet:
  # AUC of Mahalanobis distance-from-neutral, computed separately for each s>0 vs neutral
  auc_df <- d %>%
    group_by(f0, decline) %>%
    group_modify(~{
      df <- .x

      # neutral reference
      df0 <- df %>% filter(sel_num == 0)
      if (nrow(df0) < 10) {
        # not enough neutral points to define reference
        return(tibble(sel_s = character(0), auc = NA_real_))
      }

      mu <- colMeans(df0[, c("fst_hudson", "delta_tajd")], na.rm = TRUE)
      S  <- stats::cov(df0[, c("fst_hudson", "delta_tajd")], use = "complete.obs")

      # guard against singular covariance
      if (any(!is.finite(S)) || det(S) <= 0) {
        return(tibble(sel_s = character(0), auc = NA_real_))
      }

      # compute Mahalanobis distance for all points in this facet
      md_all <- stats::mahalanobis(df[, c("fst_hudson", "delta_tajd")], center = mu, cov = S)

      # for each nonzero selection level, compute AUC vs neutral
      sel_levels_here <- sort(unique(df$sel_num[df$sel_num > 0]), decreasing = FALSE)

      out <- lapply(sel_levels_here, function(sv) {
        df_sub <- df %>% filter(sel_num == 0 | sel_num == sv)

        # require enough points in both groups
        n0 <- sum(df_sub$sel_num == 0)
        n1 <- sum(df_sub$sel_num == sv)
        if (n0 < 10 || n1 < 10) {
          return(tibble(sel_s = sprintf("%.2f", sv), auc = NA_real_))
        }

        md_sub <- md_all[df$sel_num == 0 | df$sel_num == sv]
        y <- as.integer(df_sub$sel_num == sv)  # 1 = this s, 0 = neutral
        tibble(sel_s = sprintf("%.2f", sv), auc = auc_rank(md_sub, y))
      })

      bind_rows(out)
    }) %>%
    ungroup() %>%
    mutate(
      label = ifelse(is.na(auc), paste0("AUC(s=", sel_s, ")=NA"),
                     paste0("AUC(s=", sel_s, ")=", sprintf("%.2f", auc)))
    ) %>%
    group_by(f0, decline) %>%
    arrange(sel_s) %>%
    mutate(line_id = row_number()) %>%   # used to vertically stack labels
    ungroup()

  d_non0 <- d %>% filter(sel_num != 0)
  d_0    <- d %>% filter(sel_num == 0)

  # Ellipses are safe since you have lots of points, but keep a guard anyway
  ellipse_df <- d %>%
    group_by(f0, decline, sel_s) %>%
    filter(n() >= 10) %>%
    ungroup()

  # line width: neutral thicker
  ellipse_df <- ellipse_df %>%
    mutate(ell_lwd = ifelse(as.numeric(as.character(sel_s)) == 0, 0.8, 0.35))

  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +

    # ellipses (behind points)
    stat_ellipse(
      data = ellipse_df,
      aes(
        x = fst_hudson,
        y = delta_tajd,
        color = sel_s,
        group = sel_s,
        linewidth = ell_lwd
      ),
      type = "norm",
      level = 0.95,
      alpha = 0.9
    ) +
    scale_linewidth_identity() +

    # points
    geom_point(
      data = d_non0,
      aes(x = fst_hudson, y = delta_tajd, color = sel_s),
      size = 0.35, alpha = 0.8
    ) +
    geom_point(
      data = d_0,
      aes(x = fst_hudson, y = delta_tajd, color = sel_s),
      size = 0.35, alpha = 0.9
    ) +

    geom_text(
      data = auc_df,
      aes(
        x = Inf,
        y = Inf,
        label = label,
        vjust = 1.2 + 1.5 * (auc_df$line_id - 1),  # stack labels downward
      ),
      hjust = 1,
      size = 1.5,
      nudge_x = -0.01,
      inherit.aes = FALSE
    ) +


    facet_grid(
      rows = vars(f0),
      cols = vars(decline),
      labeller = labeller(
        f0 = function(x) paste0("f0=", x),
        decline = function(x) paste0("decline=", x)
      )
    ) +
    scale_color_manual(values = sel_cols, drop = FALSE) +
    theme_classic(base_size = 13) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    labs(
      x = "FST (Hudson)",
      y = expression(Delta * " Tajima's D"),
      color = "Selection (s)"
    )

  if (!is.na(xmin_pin)) p <- p + coord_cartesian(xlim = c(xmin_pin, NA))
  p
}


for (nm in names(df_split)) {
  d <- df_split[[nm]]
  if (nrow(d) == 0) next

  Ne_val  <- unique(d$Ne)
  Gen_val <- unique(d$Gen_time)

  p <- make_plot(d) +
    ggtitle(paste0("FST vs ΔTajima’s D | Ne=", Ne_val, " | Gen_time=", Gen_val))

  # output dimensions scale with facet grid size
  nrow_facets <- length(unique(d$f0))
  ncol_facets <- length(unique(d$decline))

  out_file <- file.path(out_dir, paste0("fst_vs_deltaTajd_Ne", Ne_val, "_Gen", Gen_val, ".png"))

  ggsave(
    filename = out_file,
    plot = p,
    width  = max(7, 2.4 * ncol_facets),
    height = max(6, 2.0 * nrow_facets),
    dpi = 300
  )

  message("Saved: ", out_file)
}

message("\nDone. Wrote plots to: ", out_dir)
