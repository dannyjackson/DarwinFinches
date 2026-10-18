#!/usr/bin/env Rscript
# ============================================================
# Single-species outlier power as a function of N, the number of selected
# windows in the genome. Companion to outlier_power_curves.R (which fixes
# N = 100) for the N-axis multipanel figure.
#
# Per selected condition (Ne x Gen x decline x s x f0) and each N in N_VALUES:
#   - genome of N_GENOME windows: N selected + (N_GENOME - N) neutral
#   - rank all windows by FST, take the top TOP_N (top 1%)
#   - recall = share of the N selected windows that are in the top 1%
#     (per-window detection rate; expected 0.01 with no power)
#   - P(any) = P(at least one selected window in the top 1%)
#   - repeat N_PERMUTATIONS times
#
# Run from the project directory (gsea_power_analysis):
#   Rscript ${SCRIPTDIR}/outlier_power_nsel.R
# Writes:
#   output/fst_outlier_power_nsel.tsv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
})

set.seed(123)

root_dir       <- "."
N_GENOME       <- 10000
TOP_N          <- 100
N_VALUES       <- c(1, 2, 5, 10, 20, 50, 100)
N_PERMUTATIONS <- 100
out_file       <- "output/fst_outlier_power_nsel.tsv"

read_fst <- function(file) {
  d <- read.delim(file)
  d$fst_hudson[is.finite(d$fst_hudson)]
}

all_files <- list.files(root_dir, pattern = "summary_stats\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[str_detect(all_files, "^\\./Ne_[^/]+/Gen_[0-9]+/[^/]+$")]
selected_files <- all_files[!str_detect(basename(all_files), "_s0\\.0_")]
if (length(selected_files) == 0) stop("No selected files found (run from gsea_power_analysis)")

neutral_cache <- list()
get_neutral <- function(selected_file, decline) {
  cands <- list.files(dirname(selected_file),
                      pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$", full.names = TRUE)
  if (length(cands) == 0) return(NULL)
  decl <- as.numeric(str_match(basename(cands), "^d([0-9.]+)_s")[, 2])
  hit <- which(decl == as.numeric(decline))
  if (length(hit) == 0) return(NULL)
  f <- cands[hit[1]]
  if (is.null(neutral_cache[[f]])) neutral_cache[[f]] <<- read_fst(f)
  neutral_cache[[f]]
}

one_n <- function(null_fst, sel_fst, n) {
  is_sel <- c(rep(FALSE, N_GENOME - n), rep(TRUE, n))
  res <- vapply(seq_len(N_PERMUTATIONS), function(p) {
    fst <- c(sample(null_fst, N_GENOME - n), sample(sel_fst, n))
    shuffle <- sample.int(N_GENOME)                      # random tie-breaking
    top <- shuffle[order(-fst[shuffle])][seq_len(TOP_N)]
    k <- sum(is_sel[top])
    c(recall = k / n, any = as.numeric(k > 0))
  }, numeric(2))
  c(recall = mean(res["recall", ]), prob_any = mean(res["any", ]))
}

results <- map_dfr(selected_files, function(f) {
  fn <- basename(f)
  decline <- str_match(fn, "^d([0-9.]+)_s")[, 2]
  sel_s   <- str_match(fn, "_s([0-9.]+)_af")[, 2]
  af      <- str_match(fn, "_af([0-9.]+)\\.summary")[, 2]
  Ne      <- str_match(f, "Ne_([^/]+)")[, 2]
  Gen     <- str_match(f, "Gen_([0-9]+)")[, 2]

  null_fst <- get_neutral(f, decline)
  if (is.null(null_fst) || length(null_fst) < N_GENOME - min(N_VALUES)) {
    warning("Skipping (neutral missing or too small): ", f); return(NULL)
  }
  sel_fst <- read_fst(f)
  ns <- N_VALUES[N_VALUES <= length(sel_fst) & (N_GENOME - N_VALUES) <= length(null_fst)]
  cat("Processing:", Ne, Gen, decline, sel_s, af, "\n")

  map_dfr(ns, function(n) {
    r <- one_n(null_fst, sel_fst, n)
    tibble(Ne = Ne, Gen = as.numeric(Gen), decline = as.numeric(decline),
           sel_s = as.numeric(sel_s), af = as.numeric(af), n_selected = n,
           n_permutations = N_PERMUTATIONS, recall = r[["recall"]], prob_any = r[["prob_any"]])
  })
})

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.table(results, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", out_file)
