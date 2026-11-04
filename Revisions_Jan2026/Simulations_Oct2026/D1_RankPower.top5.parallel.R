library(dplyr)
library(stringr)
library(purrr)

# ============================================================
# Parallel (three-species) version of D1_RankPower.top5.test.R,
# written to run as a Slurm array.
#
# Usage:
#   Rscript D1_RankPower.top5.parallel.R <task_id> <n_tasks>     one slice of the conditions
#   Rscript D1_RankPower.top5.parallel.R merge                   concatenate all parts
#   Rscript D1_RankPower.top5.parallel.R                         everything in one process
#
# Conditions (selected summary_stats files) are listed in a fixed
# sorted order; task k processes files with index %% n_tasks == k.
# Each condition writes its own part file to output/parallel3_parts/,
# so a killed task loses at most one condition and a rerun with the
# same arguments skips conditions already done.
#
# Same inputs as the single-species test:
#   ./Ne_*/Gen_*/d*_s*_af*.summary_stats.tsv (fst_hudson), matched to
#   the neutral file (s0.0, af0.5) of the same Ne/Gen/decline.
#
# Per permutation, N_SPECIES pseudo-species are built independently,
# each exactly as in the single-species test (N_NULL_TOTAL neutral
# windows -> N_NULL_SETS null sets of SET_SIZE plus SET_SIZE held out;
# mixed set = X selected + (SET_SIZE - X) held-out neutral windows).
# set_id is shared across species (null set k is "the same pathway"
# in every species; the mixed set has X selected windows in all
# three), while FST values are drawn independently per species.
# Each set is scored per species (mean percentile of its top-5
# windows); the cross-species score is the mean over species, the
# analogue of mean NES; the mixed set is ranked among N_SETS_TOTAL
# sets by that mean. Power = proportion of permutations in the top
# 1/5/10/25%. Also recorded: top 5% in every species separately.
#
# Merged output:
#   output/fst_rank_by_selected_windows.parallel<N_SPECIES>.nnull<N_NULL_TOTAL>.tsv
# with the single-species columns plus n_species,
# prop_all_species_top5pct and median_mean_rank_species.
# ============================================================

root_dir     <- "."
N_NULL_TOTAL <- 10000
SET_SIZE     <- 10
N_NULL_SETS  <- (N_NULL_TOTAL - SET_SIZE) / SET_SIZE
N_SETS_TOTAL <- N_NULL_SETS + 1
N_SPECIES    <- 3
X_VALUES     <- 0:SET_SIZE
N_PERM       <- 100
SEED         <- 123

stopifnot(N_NULL_SETS == round(N_NULL_SETS), N_NULL_SETS >= 1)
rank_cutoff <- function(prop) floor(prop * N_SETS_TOTAL)

part_dir     <- file.path("output", paste0("parallel", N_SPECIES, "_parts"))
results_file <- paste0("output/fst_rank_by_selected_windows.parallel", N_SPECIES,
                       ".nnull", N_NULL_TOTAL, ".tsv")

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 0) "all" else if (args[1] == "merge") "merge" else "slice"
if (mode == "slice") {
  task_id <- as.integer(args[1]); n_tasks <- as.integer(args[2])
  stopifnot(!is.na(task_id), !is.na(n_tasks), task_id >= 0, task_id < n_tasks)
}


# ------------------------------------------------------------
# One species, one permutation: per-set top-5 percentile scores
# for every X in X_VALUES (same construction as the single-species test)
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
    mixed_rank_data  <- bind_rows(selected_windows, neutral_windows) %>%
      mutate(set_type = "mixed", set_id = N_SETS_TOTAL)

    bind_rows(null_rank_data, mixed_rank_data) %>%
      mutate(fst_rank = rank(fst_hudson, ties.method = "average"),
             fst_percentile = fst_rank / n()) %>%
      group_by(set_id, set_type) %>%
      summarise(score = mean(sort(fst_percentile, decreasing = TRUE)[1:5]),
                .groups = "drop") %>%
      mutate(X_selected = X)
  })
}


rank_selected_sets_parallel <- function(null_file, selected_file, Ne, Gen, decline,
                                        sel_s, af, n_permutations = N_PERM) {

  null     <- read.delim(null_file)     %>% filter(is.finite(fst_hudson))
  selected <- read.delim(selected_file) %>% filter(is.finite(fst_hudson))

  if (nrow(null) < N_NULL_TOTAL) {
    warning(paste("Not enough usable neutral windows:", null_file, "n =", nrow(null)))
    return(NULL)
  }
  if (nrow(selected) < max(X_VALUES)) {
    warning(paste("Not enough selected windows in", selected_file))
    return(NULL)
  }

  permutation_results <- map_dfr(1:n_permutations, function(perm) {

    per_species <- map_dfr(1:N_SPECIES, function(sp) {
      score_sets_one_species(null, selected) %>% mutate(species = sp)
    }) %>%
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

    mixed_sp <- per_species %>% filter(set_type == "mixed") %>%
      group_by(X_selected) %>%
      summarise(mean_rank_species = mean(rank_species),
                max_rank_species  = max(rank_species), .groups = "drop")

    cross %>% filter(set_type == "mixed") %>%
      left_join(mixed_sp, by = "X_selected") %>%
      transmute(Ne = Ne, Gen = Gen, decline = decline, sel_s = sel_s, af = af,
                X_selected = X_selected, X_neutral = SET_SIZE - X_selected,
                permutation = perm,
                mean_percentile = mean_score,
                rank = rank_cross,
                all_species_top5 = n_species_top5 == N_SPECIES,
                mean_rank_species = mean_rank_species,
                max_rank_species = max_rank_species)
  })

  permutation_results %>%
    group_by(Ne, Gen, decline, sel_s, af, X_selected, X_neutral) %>%
    summarise(n_permutations = n(),
              median_rank = median(rank), sd_rank = sd(rank),
              min_rank = min(rank), max_rank = max(rank),
              n_sets = N_SETS_TOTAL,
              prop_top_1pct  = mean(rank <= rank_cutoff(0.01)),
              prop_top_5pct  = mean(rank <= rank_cutoff(0.05)),
              prop_top_10pct = mean(rank <= rank_cutoff(0.10)),
              prop_top_25pct = mean(rank <= rank_cutoff(0.25)),
              mean_percentile = mean(mean_percentile),
              n_species = N_SPECIES,
              prop_all_species_top5pct = mean(all_species_top5),
              median_mean_rank_species = median(mean_rank_species),
              .groups = "drop")
}


# ============================================================
# Conditions, in a fixed order shared by every task
# ============================================================
all_files <- list.files(root_dir, pattern = "summary_stats\\.tsv$", recursive = TRUE, full.names = TRUE)
all_files <- all_files[str_detect(all_files, "^\\./Ne_[^/]+/Gen_[0-9]+/[^/]+$")]
selected_files <- sort(all_files[!str_detect(basename(all_files), "_s0\\.0_")])

cond <- tibble(selected_file = selected_files) %>%
  mutate(filename = basename(selected_file),
         decline  = str_match(filename, "^d([0-9.]+)_s")[, 2],
         sel_s    = str_match(filename, "_s([0-9.]+)_af")[, 2],
         af       = str_match(filename, "_af([0-9.]+)\\.summary")[, 2],
         Ne       = str_match(selected_file, "Ne_([^/]+)")[, 2],
         Gen      = str_match(selected_file, "Gen_([0-9]+)")[, 2],
         part     = file.path(part_dir, paste0("Ne_", Ne, ".Gen_", Gen, ".", sub("\\.summary_stats\\.tsv$", "", filename), ".tsv")),
         idx      = row_number() - 1)

dir.create(part_dir, showWarnings = FALSE, recursive = TRUE)


# ============================================================
# merge: concatenate parts and stop
# ============================================================
if (mode == "merge") {
  done <- cond$part[file.exists(cond$part)]
  cat("parts found:", length(done), "of", nrow(cond), "conditions\n")
  if (length(done) < nrow(cond)) {
    cat("missing:\n"); cat(paste(" ", cond$selected_file[!file.exists(cond$part)]), sep = "\n")
  }
  all_results <- map_dfr(done, read.delim)
  write.table(all_results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", results_file, " (", nrow(all_results), " rows)")
  quit(save = "no")
}


# ============================================================
# slice / all: process this task's conditions, one part file each
# ============================================================
mine <- if (mode == "slice") cond[cond$idx %% n_tasks == task_id, ] else cond
cat("task", if (mode == "slice") task_id else "all", ":", nrow(mine), "conditions\n")

for (i in seq_len(nrow(mine))) {
  r <- mine[i, ]
  if (file.exists(r$part)) { cat("skip (done):", r$part, "\n"); next }

  neutral_candidates <- list.files(dirname(r$selected_file),
                                   pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$", full.names = TRUE)
  cand_decl <- as.numeric(str_match(basename(neutral_candidates), "^d([0-9.]+)_s")[, 2])
  hit <- which(cand_decl == as.numeric(r$decline))
  if (length(hit) == 0) { warning(paste("No matching neutral file for:", r$selected_file)); next }
  null_file <- neutral_candidates[hit[1]]

  cat("Processing:", "Ne =", r$Ne, "Gen =", r$Gen, "d =", r$decline, "s =", r$sel_s, "af =", r$af, "\n")

  # one seed per condition, derived from its position in the sorted list, so results are
  # reproducible regardless of how conditions are split across tasks
  set.seed(SEED + r$idx)

  res <- rank_selected_sets_parallel(
    null_file = null_file, selected_file = r$selected_file,
    Ne = r$Ne, Gen = as.numeric(r$Gen), decline = as.numeric(r$decline),
    sel_s = as.numeric(r$sel_s), af = as.numeric(r$af), n_permutations = N_PERM)
  if (is.null(res)) next

  tmp <- paste0(r$part, ".tmp")
  write.table(res, tmp, sep = "\t", quote = FALSE, row.names = FALSE)
  file.rename(tmp, r$part)                 # a part file is only ever complete
}

if (mode == "all") {
  all_results <- map_dfr(cond$part[file.exists(cond$part)], read.delim)
  write.table(all_results, results_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("wrote ", results_file)
}
message("Done.")
