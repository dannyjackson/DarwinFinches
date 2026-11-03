library(dplyr)
library(stringr)
library(purrr)

set.seed(123)

root_dir <- "."

# ============================================================
# Parallel (three-species) version of D1_RankPower.top5.test.R
#
# Same inputs: ./Ne_*/Gen_*/d*_s*_af*.summary_stats.tsv (fst_hudson),
# matched to the neutral file (s0.0, af0.5) of the same Ne/Gen/decline.
#
# Per permutation, N_SPECIES pseudo-species are built independently,
# each exactly as in the single-species test:
#   N_NULL_TOTAL - SET_SIZE neutral windows -> N_NULL_SETS null sets
#   SET_SIZE held-out neutral windows        -> fill the mixed set
#   mixed set = X selected windows + (SET_SIZE - X) held-out neutral
# Set identities (set_id) are shared across species, like gene
# identities in the parallel outlier test: null set k in species 1
# is "the same pathway" as null set k in species 2 and 3, but its
# windows are drawn independently in each species. The mixed set
# has X selected windows in every species (parallel selection), with
# FST values drawn independently per species from the same selected
# replicates.
#
# Each set is scored per species as in the single-species test
# (mean percentile of its top-5 windows); the cross-species score is
# the mean of the per-species scores, the analogue of mean NES.
# Power = proportion of permutations in which the mixed set ranks in
# the top 1/5/10/25% of N_SETS_TOTAL sets by cross-species score.
# Also recorded: the mixed set's per-species rank summary and the
# proportion of permutations in which it is in the top 5% in every
# species (the "outlier in all species" criterion applied to sets).
#
# Output: output/fst_rank_by_selected_windows.parallel<N_SPECIES>.nnull<N_NULL_TOTAL>.tsv
# Columns match the single-species output, plus the extra columns
# listed at the end of rank_selected_sets_parallel().
# ============================================================

N_NULL_TOTAL <- 10000
SET_SIZE     <- 10
N_NULL_SETS  <- (N_NULL_TOTAL - SET_SIZE) / SET_SIZE
N_SETS_TOTAL <- N_NULL_SETS + 1   # null sets + the mixed set
N_SPECIES    <- 3
X_VALUES     <- 0:SET_SIZE        # selected windows in the mixed set (same in every species)

stopifnot(N_NULL_SETS == round(N_NULL_SETS), N_NULL_SETS >= 1)

rank_cutoff <- function(prop) floor(prop * N_SETS_TOTAL)


# ------------------------------------------------------------
# One species, one permutation: per-set top-5 percentile scores
# for every X in X_VALUES. Returns a tibble: X_selected, set_id,
# set_type, score. Identical construction to the single-species
# test (shuffle, split, mixed set, rank, mean of top 5).
# ------------------------------------------------------------
score_sets_one_species <- function(null, selected) {

  null_shuffled   <- null %>% slice_sample(n = N_NULL_TOTAL)
  null_background <- null_shuffled[1:(N_NULL_TOTAL - SET_SIZE), ] %>%
    mutate(null_group = rep(1:N_NULL_SETS, each = SET_SIZE))
  null_holdout    <- null_shuffled[(N_NULL_TOTAL - SET_SIZE + 1):N_NULL_TOTAL, ]

  null_rank_data <- null_background %>%
    mutate(set_type = "null", set_id = null_group)

  map_dfr(X_VALUES, function(X) {

    selected_windows <- selected %>% slice_sample(n = X)
    neutral_windows  <- null_holdout %>% slice_sample(n = SET_SIZE - X)

    mixed_rank_data <- bind_rows(selected_windows, neutral_windows) %>%
      mutate(set_type = "mixed", set_id = N_SETS_TOTAL)

    bind_rows(null_rank_data, mixed_rank_data) %>%
      mutate(
        fst_rank = rank(fst_hudson, ties.method = "average"),
        fst_percentile = fst_rank / n()
      ) %>%
      group_by(set_id, set_type) %>%
      summarise(
        score = mean(sort(fst_percentile, decreasing = TRUE)[1:5]),
        .groups = "drop"
      ) %>%
      mutate(X_selected = X)
  })
}


rank_selected_sets_parallel <- function(null_file, selected_file, Ne, Gen, decline,
                                        sel_s, af, n_permutations = 100) {

  null <- read.delim(null_file) %>%
    filter(is.finite(fst_hudson))

  selected <- read.delim(selected_file) %>%
    filter(is.finite(fst_hudson))

  if (nrow(null) < N_NULL_TOTAL) {
    warning(paste("Not enough usable neutral windows:",
                  null_file, "n =", nrow(null)))
    return(NULL)
  }
  if (nrow(selected) < max(X_VALUES)) {
    warning(paste("Not enough selected windows in", selected_file))
    return(NULL)
  }

  permutation_results <- map_dfr(1:n_permutations, function(perm) {

    # Independent construction in each species; set_id shared
    per_species <- map_dfr(1:N_SPECIES, function(sp) {
      score_sets_one_species(null, selected) %>% mutate(species = sp)
    })

    # Per-species rank of every set (1 = strongest), then cross-species mean score
    per_species <- per_species %>%
      group_by(X_selected, species) %>%
      mutate(rank_species = rank(-score, ties.method = "min")) %>%
      ungroup()

    cross <- per_species %>%
      group_by(X_selected, set_id, set_type) %>%
      summarise(mean_score = mean(score),
                n_species_top5 = sum(rank_species <= rank_cutoff(0.05)),
                .groups = "drop") %>%
      group_by(X_selected) %>%
      mutate(rank_cross = rank(-mean_score, ties.method = "min")) %>%
      ungroup()

    mixed_cross <- cross %>% filter(set_type == "mixed")
    mixed_sp    <- per_species %>% filter(set_type == "mixed") %>%
      group_by(X_selected) %>%
      summarise(mean_rank_species = mean(rank_species),
                max_rank_species  = max(rank_species),
                .groups = "drop")

    mixed_cross %>%
      left_join(mixed_sp, by = "X_selected") %>%
      transmute(
        Ne = Ne, Gen = Gen, decline = decline, sel_s = sel_s, af = af,
        X_selected = X_selected,
        X_neutral  = SET_SIZE - X_selected,
        permutation = perm,
        mean_percentile = mean_score,       # cross-species mean of top-5 percentile scores
        rank = rank_cross,                  # rank among N_SETS_TOTAL sets by cross-species score
        all_species_top5 = n_species_top5 == N_SPECIES,
        mean_rank_species = mean_rank_species,
        max_rank_species = max_rank_species
      )
  })

  permutation_results %>%
    group_by(Ne, Gen, decline, sel_s, af, X_selected, X_neutral) %>%
    summarise(
      n_permutations = n(),
      median_rank = median(rank),
      sd_rank = sd(rank),
      min_rank = min(rank),
      max_rank = max(rank),
      n_sets = N_SETS_TOTAL,
      prop_top_1pct  = mean(rank <= rank_cutoff(0.01)),
      prop_top_5pct  = mean(rank <= rank_cutoff(0.05)),
      prop_top_10pct = mean(rank <= rank_cutoff(0.10)),
      prop_top_25pct = mean(rank <= rank_cutoff(0.25)),
      mean_percentile = mean(mean_percentile),
      # extra columns (not in the single-species output)
      n_species = N_SPECIES,
      prop_all_species_top5pct = mean(all_species_top5),   # top 5% in every species separately
      median_mean_rank_species = median(mean_rank_species), # typical single-species rank of the mixed set
      .groups = "drop"
    )
}


# ============================================================
# Find all summary-statistic files (as in the single-species test)
# ============================================================

all_files <- list.files(root_dir, pattern = "summary_stats\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[str_detect(all_files, "^\\./Ne_[^/]+/Gen_[0-9]+/[^/]+$")]
selected_files <- all_files[!str_detect(basename(all_files), "_s0\\.0_")]


all_results <- map_dfr(selected_files, function(selected_file) {

  filename <- basename(selected_file)
  decline <- str_match(filename, "^d([0-9.]+)_s")[, 2]
  sel_s   <- str_match(filename, "_s([0-9.]+)_af")[, 2]
  af      <- str_match(filename, "_af([0-9.]+)\\.summary")[, 2]
  Ne      <- str_match(selected_file, "Ne_([^/]+)")[, 2]
  Gen     <- str_match(selected_file, "Gen_([0-9]+)")[, 2]
  selected_dir <- dirname(selected_file)

  neutral_candidates <- list.files(selected_dir,
                                   pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$",
                                   full.names = TRUE)
  if (length(neutral_candidates) == 0) {
    warning(paste("No neutral file found for:", selected_file)); return(NULL)
  }
  candidate_declines <- as.numeric(str_match(basename(neutral_candidates), "^d([0-9.]+)_s")[, 2])
  match_index <- which(candidate_declines == as.numeric(decline))
  if (length(match_index) == 0) {
    warning(paste("No matching neutral decline file for:", selected_file)); return(NULL)
  }
  null_file <- neutral_candidates[match_index[1]]

  cat("Processing:", "Ne =", Ne, "Gen =", Gen, "d =", decline,
      "s =", sel_s, "af =", af, "\n")

  rank_selected_sets_parallel(
    null_file = null_file, selected_file = selected_file,
    Ne = Ne, Gen = as.numeric(Gen), decline = as.numeric(decline),
    sel_s = as.numeric(sel_s), af = as.numeric(af),
    n_permutations = 100
  )
})


# ============================================================
# Save results
# ============================================================

dir.create("output", showWarnings = FALSE)

results_file <- paste0("output/fst_rank_by_selected_windows.parallel", N_SPECIES,
                       ".nnull", N_NULL_TOTAL, ".tsv")

write.table(all_results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)

print(all_results)
message("\nDone. Results written to: ", results_file)
