#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

# ============================================================
# Settings
# ============================================================

root_dir <- "/xdisk/mcnew/finches/dannyjackson/simulations/power_analysis"

out_dir <- file.path(
  root_dir,
  "plots_mean_comparisons_neutral_vs_selected"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Significance threshold
alpha <- 0.05

# Which p-value determines the binary heatmap?
#
# "p_value"      = raw Welch t-test p-value
# "p_adj_stat"   = BH correction separately within each statistic
# "p_adj_global" = BH correction across all tests
significance_column <- "p_adj_stat"

# Optional filters
keep_ne      <- NULL
keep_gen     <- NULL
keep_decline <- NULL
keep_f0      <- NULL
keep_sel     <- NULL


# ============================================================
# Discover files
# ============================================================

files <- Sys.glob(
  file.path(root_dir, "Ne_*", "Gen_*", "*.summary_stats.tsv")
)

message("Found ", length(files), " summary_stats.tsv files")

if (length(files) == 0) {
  stop("No *.summary_stats.tsv files found under: ", root_dir)
}


# ============================================================
# Parse metadata
# ============================================================

parse_ne <- function(ne_chr) {

  x <- str_remove(ne_chr, "^Ne_")

  if (str_detect(x, "k$")) {
    as.numeric(str_remove(x, "k$")) * 1000
  } else {
    suppressWarnings(as.numeric(x))
  }
}


parse_from_filename <- function(fname) {

  base <- basename(fname)

  # Expected:
  # d0.80_s0.05_af0.05.summary_stats.tsv
  #
  # Also accepts:
  # d1.00_s0.0_af0.5.summary_stats.tsv

  m <- str_match(
    base,
    "^d([0-9]+(?:\\.[0-9]+)?)_s([0-9]+(?:\\.[0-9]+)?)_af([0-9]+(?:\\.[0-9]+)?)\\.summary_stats\\.tsv$"
  )

  if (is.na(m[1, 1])) {
    return(
      list(
        decline = NA_real_,
        sel_s   = NA_real_,
        f0      = NA_real_
      )
    )
  }

  list(
    decline = as.numeric(m[1, 2]),
    sel_s   = as.numeric(m[1, 3]),
    f0      = as.numeric(m[1, 4])
  )
}


parse_from_path <- function(path) {

  rel <- gsub(
    paste0("^", root_dir, "/"),
    "",
    path
  )

  parts <- str_split(rel, "/", simplify = TRUE)

  # expected:
  # Ne_100k/Gen_15/d0.80_s0.05_af0.5.summary_stats.tsv

  if (ncol(parts) < 3) {
    return(tibble())
  }

  ne_chr  <- parts[1]
  gen_chr <- parts[2]
  fn      <- parts[3]

  fp <- parse_from_filename(fn)

  tibble(
    file     = path,
    Ne       = parse_ne(ne_chr),
    Gen_time = suppressWarnings(
      as.numeric(str_remove(gen_chr, "^Gen_"))
    ),
    decline  = fp$decline,
    f0       = fp$f0,
    sel_s    = fp$sel_s
  )
}


meta <- map_dfr(files, parse_from_path) %>%
  filter(
    !is.na(Ne),
    !is.na(Gen_time),
    !is.na(decline),
    !is.na(sel_s)
  )


if (nrow(meta) == 0) {
  stop("Parsed 0 metadata rows.")
}


message("Parsed metadata for ", nrow(meta), " files")

message("Unique parameter values:")
message("  Ne:        ", paste(sort(unique(meta$Ne)), collapse = ", "))
message("  Gen_time:  ", paste(sort(unique(meta$Gen_time)), collapse = ", "))
message("  decline:   ", paste(sort(unique(meta$decline)), collapse = ", "))
message("  f0:        ", paste(sort(unique(meta$f0)), collapse = ", "))
message("  sel_s:     ", paste(sort(unique(meta$sel_s)), collapse = ", "))


# ============================================================
# Optional filtering
# ============================================================

if (!is.null(keep_ne)) {
  meta <- meta %>% filter(Ne %in% keep_ne)
}

if (!is.null(keep_gen)) {
  meta <- meta %>% filter(Gen_time %in% keep_gen)
}

if (!is.null(keep_decline)) {
  meta <- meta %>% filter(decline %in% keep_decline)
}

if (!is.null(keep_sel)) {
  # Always retain neutral files.
  meta <- meta %>%
    filter(sel_s == 0 | sel_s %in% keep_sel)
}

if (!is.null(keep_f0)) {
  # f0 matters for selected runs.
  # Keep neutral runs regardless of their stored f0 label.
  meta <- meta %>%
    filter(sel_s == 0 | f0 %in% keep_f0)
}

if (nrow(meta) == 0) {
  stop("After filtering, no files remain.")
}


# ============================================================
# Read summary-statistic files
# ============================================================

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


read_one <- function(
    file,
    Ne,
    Gen_time,
    decline,
    f0,
    sel_s
) {

  lines <- readLines(file, warn = FALSE)

  if (length(lines) < 2) {
    return(tibble())
  }

  expected_n <- length(
    strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  )

  good <- vapply(
    strsplit(lines, "\t", fixed = TRUE),
    length,
    integer(1)
  ) == expected_n

  df <- readr::read_tsv(
    I(lines[good]),
    col_types = forced_types,
    show_col_types = FALSE
  )

  need <- c(
    "fst_hudson",
    "delta_tajd"
  )

  missing <- setdiff(need, names(df))

  if (length(missing) > 0) {
    stop(
      "Missing required columns in ",
      file,
      ": ",
      paste(missing, collapse = ", ")
    )
  }

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

message(
  "Finished reading data. Total rows loaded: ",
  nrow(df_all)
)

if (nrow(df_all) == 0) {
  stop("No rows read from files.")
}


df_all <- df_all %>%
  mutate(
    sel_num = as.numeric(sel_s)
  )


# ============================================================
# Welch t-test helper
# ============================================================

safe_ttest <- function(neutral, selected) {

  neutral <- neutral[
    is.finite(neutral)
  ]

  selected <- selected[
    is.finite(selected)
  ]

  n_neutral  <- length(neutral)
  n_selected <- length(selected)

  mean_neutral <- if (n_neutral > 0) {
    mean(neutral)
  } else {
    NA_real_
  }

  mean_selected <- if (n_selected > 0) {
    mean(selected)
  } else {
    NA_real_
  }

  difference <- mean_selected - mean_neutral

  if (n_neutral < 2 || n_selected < 2) {

    return(
      tibble(
        n_neutral     = n_neutral,
        n_selected    = n_selected,
        mean_neutral  = mean_neutral,
        mean_selected = mean_selected,
        difference    = difference,
        t_statistic   = NA_real_,
        df            = NA_real_,
        p_value       = NA_real_
      )
    )
  }

  test <- tryCatch(
    t.test(
      x = selected,
      y = neutral,
      alternative = "two.sided",
      var.equal = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(test)) {

    return(
      tibble(
        n_neutral     = n_neutral,
        n_selected    = n_selected,
        mean_neutral  = mean_neutral,
        mean_selected = mean_selected,
        difference    = difference,
        t_statistic   = NA_real_,
        df            = NA_real_,
        p_value       = NA_real_
      )
    )
  }

  tibble(
    n_neutral     = n_neutral,
    n_selected    = n_selected,
    mean_neutral  = mean_neutral,
    mean_selected = mean_selected,
    difference    = difference,
    t_statistic   = unname(test$statistic),
    df            = unname(test$parameter),
    p_value       = test$p.value
  )
}


# ============================================================
# Neutral-vs-selected tests
#
# IMPORTANT:
# This follows the original power-analysis logic.
#
# Neutral simulations are pooled within:
#   Ne x Gen_time x decline
#
# Their stored f0 label is ignored because f0 is irrelevant
# when s = 0.
#
# Each selected:
#   s x f0
#
# combination is then compared to that neutral pool.
#
# Composite statistic matches the original analysis:
#
#   composite = z(FST) - z(delta Tajima's D)
#
# where the z-scores are computed within the combined
# neutral + selected comparison.
# ============================================================

group_vars <- c(
  "Ne",
  "Gen_time",
  "decline"
)


test_results <- df_all %>%

  group_by(
    across(all_of(group_vars))
  ) %>%

  group_modify(~{

    d <- .x

    # --------------------------------
    # Neutral reference simulations
    # --------------------------------

    neutral_pool <- d %>%
      filter(sel_num == 0)

    if (nrow(neutral_pool) == 0) {
      return(tibble())
    }

    # --------------------------------
    # Selected simulations
    # --------------------------------

    selected_dat <- d %>%
      filter(
        sel_num > 0,
        !is.na(f0)
      )

    if (nrow(selected_dat) == 0) {
      return(tibble())
    }

    # Every selected s x f0 combination
    combinations <- selected_dat %>%
      distinct(
        sel_num,
        f0
      ) %>%
      arrange(
        f0,
        sel_num
      )

    pmap_dfr(
      combinations,
      function(sel_num, f0) {

        selected_block <- selected_dat %>%
          filter(
            .data$sel_num == sel_num,
            .data$f0 == f0
          )

        # --------------------------------
        # 1. Raw FST
        # --------------------------------

        fst_test <- safe_ttest(
          neutral_pool$fst_hudson,
          selected_block$fst_hudson
        ) %>%
          mutate(
            statistic = "fst_hudson",
            statistic_label = "FST (Hudson)",
            .before = 1
          )

        # --------------------------------
        # 2. Raw delta Tajima's D
        # --------------------------------

        taj_test <- safe_ttest(
          neutral_pool$delta_tajd,
          selected_block$delta_tajd
        ) %>%
          mutate(
            statistic = "delta_tajd",
            statistic_label = "Delta Tajima's D",
            .before = 1
          )

        # --------------------------------
        # 3. Composite statistic
        #
        # Match original power analysis:
        # z-scores are calculated within
        # neutral + this one selected block.
        # --------------------------------

        dd <- bind_rows(
          neutral_pool %>%
            mutate(test_group = "neutral"),
          selected_block %>%
            mutate(test_group = "selected")
        )

        dd <- dd %>%
          mutate(
            z_fst = as.numeric(
              scale(fst_hudson)
            ),
            z_taj = as.numeric(
              scale(delta_tajd)
            ),
            composite = z_fst - z_taj
          )

        composite_test <- safe_ttest(
          dd$composite[
            dd$test_group == "neutral"
          ],
          dd$composite[
            dd$test_group == "selected"
          ]
        ) %>%
          mutate(
            statistic = "composite",
            statistic_label = "Composite: z(FST) - z(Delta Tajima's D)",
            .before = 1
          )

        bind_rows(
          fst_test,
          taj_test,
          composite_test
        ) %>%
          mutate(
            sel_s = sel_num,
            f0_selected = f0,
            .before = 1
          )
      }
    )
  }) %>%

  ungroup()


if (nrow(test_results) == 0) {
  stop(
    "No neutral-vs-selected tests were produced. ",
    "Check that neutral and selected simulations are both present."
  )
}


# ============================================================
# Multiple-testing correction
# ============================================================

# Correct separately within each statistic.
test_results <- test_results %>%

  group_by(statistic) %>%

  mutate(
    p_adj_stat = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%

  ungroup()


# Also provide one correction across every test.
test_results <- test_results %>%

  mutate(
    p_adj_global = p.adjust(
      p_value,
      method = "BH"
    )
  )


# ============================================================
# Define binary significance
# ============================================================

if (!significance_column %in% names(test_results)) {
  stop(
    "significance_column must be one of: ",
    paste(
      c(
        "p_value",
        "p_adj_stat",
        "p_adj_global"
      ),
      collapse = ", "
    )
  )
}


test_results <- test_results %>%

  mutate(
    p_for_plot = .data[[significance_column]],

    significant =
      !is.na(p_for_plot) &
      p_for_plot < alpha
  )


# ============================================================
# Save full statistical table
# ============================================================

out_table <- file.path(
  out_dir,
  "neutral_vs_selected_ttests_with_composite.tsv"
)

write_tsv(
  test_results,
  out_table
)

message(
  "Wrote test results: ",
  out_table
)


# ============================================================
# Print summary
# ============================================================

test_results %>%

  group_by(
    statistic,
    statistic_label
  ) %>%

  summarise(
    n_tests = n(),
    n_significant = sum(
      significant,
      na.rm = TRUE
    ),
    proportion_significant = mean(
      significant,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%

  print()


# ============================================================
# Prepare plotting variables
# ============================================================

plot_df <- test_results %>%

  mutate(

    Ne = factor(
      Ne,
      levels = sort(
        unique(Ne)
      )
    ),

    decline = factor(
      decline,
      levels = sort(
        unique(decline)
      )
    ),

    Gen_time = factor(
      Gen_time,
      levels = sort(
        unique(Gen_time)
      )
    ),

    sel_s = factor(
      sel_s,
      levels = sort(
        unique(sel_s)
      )
    ),

    f0_selected = factor(
      f0_selected,
      levels = sort(
        unique(f0_selected)
      )
    )
  )


# ============================================================
# Binary heatmaps
#
# x     = selection coefficient
# y     = generations between samples
# rows  = decline
# cols  = Ne
#
# Separate plot for each:
#   statistic x f0
# ============================================================

combos <- plot_df %>%

  distinct(
    f0_selected,
    statistic,
    statistic_label
  )


for (i in seq_len(nrow(combos))) {

  f0_val <- combos$f0_selected[i]

  stat_val <- combos$statistic[i]

  stat_label <- combos$statistic_label[i]

  d <- plot_df %>%

    filter(
      f0_selected == f0_val,
      statistic == stat_val
    )

  p <- ggplot(
    d,
    aes(
      x = sel_s,
      y = Gen_time,
      fill = significant
    )
  ) +

    geom_tile(
      color = "white",
      linewidth = 0.4
    ) +

    facet_grid(
      decline ~ Ne,
      labeller = label_both
    ) +

    scale_fill_manual(
      values = c(
        "FALSE" = "grey90",
        "TRUE"  = "firebrick3"
      ),
      labels = c(
        "FALSE" = "Not significant",
        "TRUE"  = "Significant"
      ),
      name = NULL,
      na.value = "grey50"
    ) +

    labs(
      x = "Selection coefficient (s)",
      y = "Generations between samples",

      title = paste0(
        stat_label,
        ": selected vs neutral | f0 = ",
        f0_val
      ),

      subtitle = paste0(
        "Welch t-test; significance based on ",
        significance_column,
        " < ",
        alpha
      )
    ) +

    theme_bw(
      base_size = 13
    ) +

    theme(
      panel.grid = element_blank(),
      strip.background =
        element_rect(fill = "grey90"),
      legend.position = "bottom"
    )


  outfile <- file.path(
    out_dir,
    paste0(
      "ttest_significance.",
      stat_val,
      ".f0_",
      f0_val,
      ".png"
    )
  )


  ggsave(
    outfile,
    plot = p,
    width = 10,
    height = 3,
    dpi = 300
  )


  message(
    "Saved: ",
    outfile
  )
}


message(
  "\nDone. Results written to: ",
  out_dir
)