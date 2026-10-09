library(dplyr)
library(stringr)
library(purrr)
library(ggplot2)

set.seed(123)

root_dir <- "."

# ============================================================
# Function: rank mixed selected/null gene sets
#
# For each X = 1:9:
#   - repeat the shuffle/ranking procedure n_permutations times
#   - rank all windows by FST
#   - summarize each 10-window set by mean percentile rank
#   - rank the mixed set among 999 neutral sets + itself
#
# Within each permutation:
#   - 9990 neutral windows -> 999 background sets
#   - 10 neutral windows -> held out for mixed set
# ============================================================


rank_selected_sets <- function(null_file, selected_file, Ne, Gen, decline,
                               sel_s, af, n_permutations = 100) {

  null <- read.delim(null_file) %>%
    filter(is.finite(fst_hudson))

  selected <- read.delim(selected_file) %>%
    filter(is.finite(fst_hudson))

  # Need 9990 neutral windows for background + 10 held out
  if (nrow(null) < 10000) {
    warning(paste("Not enough usable neutral windows:",
                  null_file, "n =", nrow(null)))
    return(NULL)
  }

  permutation_results <- map_dfr(1:n_permutations, function(perm) {

    # Shuffle and split neutral windows
    null_shuffled <- null %>% slice_sample(n = 10000)
    null_background <- null_shuffled[1:9990, ]
    null_holdout <- null_shuffled[9991:10000, ]

    # Create 999 non-overlapping neutral groups of 10
    null_background <- null_background %>%
      mutate(null_group = rep(1:999, each = 10))

    # Test X = 0:10 selected windows
    map_dfr(0:10, function(X) {

      if (nrow(selected) < X) {
        warning(paste("Not enough selected windows for X =", X,
                      "in", selected_file))
        return(NULL)
      }

      # Build mixed set
      selected_windows <- selected %>% slice_sample(n = X)
      neutral_windows <- null_holdout %>% slice_sample(n = 10 - X)

      mixed_windows <- bind_rows(selected_windows, neutral_windows) %>%
        mutate(set_type = "mixed")

      # Label null and mixed sets
      null_rank_data <- null_background %>%
        mutate(set_type = "null", set_id = null_group)

      mixed_rank_data <- mixed_windows %>%
        mutate(set_id = 1000)

      # Rank all 10,000 windows together by FST
      ranked_windows <- bind_rows(null_rank_data, mixed_rank_data) %>%
        mutate(
          fst_rank = rank(fst_hudson, ties.method = "average"),
          fst_percentile = fst_rank / n()
        )

      # Mean percentile rank for each 10-window set
      set_scores <- ranked_windows %>%
        group_by(set_id, set_type) %>%
        summarise(
          top5_percentile = mean(sort(fst_percentile, decreasing = TRUE)[1:5]),
          .groups = "drop"
        )

      # Score and rank of focal mixed set
      mixed_score <- set_scores %>%
        filter(set_type == "mixed") %>%
        pull(top5_percentile)

      mixed_rank <- rank(
        -set_scores$top5_percentile,
        ties.method = "min"
      )[set_scores$set_type == "mixed"]

      tibble(
        Ne = Ne,
        Gen = Gen,
        decline = decline,
        sel_s = sel_s,
        af = af,
        X_selected = X,
        X_neutral = 10 - X,
        permutation = perm,
        mean_percentile = mixed_score,
        rank = mixed_rank
      )
    })
  })

  # Summarize across permutations
  results <- permutation_results %>%
    group_by(Ne, Gen, decline, sel_s, af, X_selected, X_neutral) %>%
    summarise(
      n_permutations = n(),
      median_rank = median(rank),
      sd_rank = sd(rank),
      min_rank = min(rank),
      max_rank = max(rank),
      prop_top_1pct = mean(rank <= 10),
      prop_top_5pct = mean(rank <= 50),
      prop_top_10pct = mean(rank <= 100),
      prop_top_25pct = mean(rank <= 250),
      mean_percentile = mean(mean_percentile),
      .groups = "drop"
    )

  return(results)
}



# ============================================================
# Find all summary-statistic files
# ============================================================

all_files <- list.files(
  root_dir,
  pattern = "summary_stats\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Keep only positively selected files
selected_files <- all_files[
  !str_detect(
    basename(all_files),
    "_s0\\.0_"
  )
]


# ============================================================
# Process each selected file
# ============================================================

all_results <- map_dfr(selected_files, function(selected_file) {

  filename <- basename(selected_file)

  # Extract parameters from filename
  # Example:
  # d0.80_s0.05_af0.05.summary_stats.tsv

  decline <- str_match(
    filename,
    "^d([0-9.]+)_s"
  )[, 2]

  sel_s <- str_match(
    filename,
    "_s([0-9.]+)_af"
  )[, 2]

  af <- str_match(
    filename,
    "_af([0-9.]+)\\.summary"
  )[, 2]

  # Extract Ne and generation time from path
  # Example:
  # ./Ne_100k/Gen_5/...

  Ne <- str_match(
    selected_file,
    "Ne_([^/]+)"
  )[, 2]

  Gen <- str_match(
    selected_file,
    "Gen_([0-9]+)"
  )[, 2]

  selected_dir <- dirname(selected_file)

  # ----------------------------------------------------------
  # Find matching neutral file
  #
  # Neutral files have:
  #   s0.0
  #   af0.5
  # ----------------------------------------------------------

  neutral_candidates <- list.files(
    selected_dir,
    pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$",
    full.names = TRUE
  )

  if (length(neutral_candidates) == 0) {
    warning(paste("No neutral file found for:", selected_file))
    return(NULL)
  }

  # Compare decline numerically so 0.8 == 0.80
  candidate_declines <- as.numeric(
    str_match(
      basename(neutral_candidates),
      "^d([0-9.]+)_s"
    )[, 2]
  )

  match_index <- which(
    candidate_declines == as.numeric(decline)
  )

  if (length(match_index) == 0) {
    warning(paste(
      "No matching neutral decline file for:",
      selected_file
    ))
    return(NULL)
  }

  null_file <- neutral_candidates[match_index[1]]

  cat(
    "Processing:",
    "Ne =", Ne,
    "Gen =", Gen,
    "d =", decline,
    "s =", sel_s,
    "af =", af,
    "\n"
  )

  rank_selected_sets(
    null_file = null_file,
    selected_file = selected_file,
    Ne = Ne,
    Gen = as.numeric(Gen),
    decline = as.numeric(decline),
    sel_s = as.numeric(sel_s),
    af = as.numeric(af),
    n_permutations = 100
  )
})


# ============================================================
# Save results
# ============================================================

dir.create(
  "output",
  showWarnings = FALSE
)

write.table(
  all_results,
  "output/fst_rank_by_selected_windows.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

print(all_results)


# ============================================================
# Prepare plotting variables
# ============================================================

plot_df <- all_results %>%
  mutate(

    # Order Ne numerically despite labels such as 100k
    Ne_num = case_when(
      str_detect(Ne, "k$") ~
        as.numeric(str_remove(Ne, "k$")) * 1000,
      TRUE ~ as.numeric(Ne)
    ),

    Ne = factor(
      Ne,
      levels = unique(
        Ne[order(Ne_num)]
      )
    ),

    decline = factor(
      decline,
      levels = sort(unique(decline))
    ),

    Gen = factor(
      Gen,
      levels = sort(unique(Gen))
    ),

    sel_s = factor(
      sel_s,
      levels = sort(unique(sel_s))
    ),

    af = factor(
      af,
      levels = sort(unique(af))
    ),

    X_selected = factor(
      X_selected,
      levels = 1:9
    )
  )


# ============================================================
# Plotting structure
#
# x    = selection coefficient
# y    = number of selected windows
# rows = decline
# cols = Ne
#
# Separate plots for each:
#   starting allele frequency x generation time
# ============================================================

plot_dir <- file.path(
  "output",
  "fst_rank_heatmaps"
)

dir.create(
  plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

combos <- plot_df %>%
  distinct(
    af,
    Gen
  ) %>%
  arrange(
    af,
    Gen
  )


# ============================================================
# 1. Median-rank heatmaps
#
# Rank 1 = strongest enrichment
# Rank 1000 = weakest enrichment
# ============================================================

for (i in seq_len(nrow(combos))) {

  af_val <- combos$af[i]
  gen_val <- combos$Gen[i]

  d <- plot_df %>%
    filter(
      af == af_val,
      Gen == gen_val
    )

  p <- ggplot(
    d,
    aes(
      x = sel_s,
      y = X_selected,
      fill = median_rank
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
    scale_fill_gradient(
      low = "firebrick3",
      high = "grey90",
      limits = c(1, 1000),
      name = "Median rank"
    ) +
    labs(
      x = "Selection coefficient (s)",
      y = "Number of selected windows",
      title = paste0(
        "Median pathway rank | f0 = ",
        af_val,
        ", generations = ",
        gen_val
      ),
      subtitle = paste0(
        "Median rank across permutations; ",
        "rank 1 = strongest enrichment"
      )
    ) +
    theme_bw(
      base_size = 13
    ) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "bottom"
    )

  outfile <- file.path(
    plot_dir,
    paste0(
      "median_rank.",
      "af_",
      af_val,
      ".Gen_",
      gen_val,
      ".png"
    )
  )

  ggsave(
    outfile,
    plot = p,
    width = 10,
    height = 7,
    dpi = 300
  )

  message("Saved: ", outfile)
}


# ============================================================
# 2. Proportion-in-top-X heatmaps
#
# Each plot shows the proportion of permutations in which
# the mixed pathway reached a specified part of the ranking.
#
# Four separate metrics:
#   prop_top_1pct
#   prop_top_5pct
#   prop_top_10pct
#   prop_top_25pct
# ============================================================

rank_metrics <- tibble(
  variable = c(
    "prop_top_1pct",
    "prop_top_5pct",
    "prop_top_10pct",
    "prop_top_25pct"
  ),
  label = c(
    "Top 1%",
    "Top 5%",
    "Top 10%",
    "Top 25%"
  )
)


for (m in seq_len(nrow(rank_metrics))) {

  metric <- rank_metrics$variable[m]
  metric_label <- rank_metrics$label[m]

  for (i in seq_len(nrow(combos))) {

    af_val <- combos$af[i]
    gen_val <- combos$Gen[i]

    d <- plot_df %>%
      filter(
        af == af_val,
        Gen == gen_val
      )

    p <- ggplot(
      d,
      aes(
        x = sel_s,
        y = X_selected,
        fill = .data[[metric]]
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
      scale_fill_gradient(
        low = "white",
        high = "firebrick3",
        limits = c(0, 1),
        name = "Proportion"
      ) +
      labs(
        x = "Selection coefficient (s)",
        y = "Number of selected windows",
        title = paste0(
          "Probability of ranking in ",
          metric_label,
          " | f0 = ",
          af_val,
          ", generations = ",
          gen_val
        ),
        subtitle = paste0(
          "Proportion of permutations in which the mixed ",
          "10-window set reached the ",
          metric_label,
          " of pathway ranks"
        )
      ) +
      theme_bw(
        base_size = 13
      ) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        legend.position = "bottom"
      )

    outfile <- file.path(
      plot_dir,
      paste0(
        metric,
        ".af_",
        af_val,
        ".Gen_",
        gen_val,
        ".png"
      )
    )

    ggsave(
      outfile,
      plot = p,
      width = 10,
      height = 7,
      dpi = 300
    )

    message("Saved: ", outfile)
  }
}


message(
  "\nDone. Results and heatmaps written to: ",
  plot_dir
)