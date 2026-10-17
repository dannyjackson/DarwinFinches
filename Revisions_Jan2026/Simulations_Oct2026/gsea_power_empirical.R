#!/usr/bin/env Rscript
# ============================================================
# Power of the empirical GSEA pipeline to detect parallel selection
# on one GO term (default: blood coagulation, GO:0007596)
#
# For ONE simulated condition (Ne, Gen, decline, s, f0), repeat N_ITER times:
#   For each species (CRA, FOR, PAR), independently:
#     - start from that species' empirical gene list
#       (gene_scores/<sp>.fst.gene_scores.tsv)
#     - give every gene an FST drawn from the matching NEUTRAL simulations
#       (same Ne, Gen, decline)
#     - give N genes in the focal GO term an FST drawn from the SELECTED
#       simulations, for every N in n_sel (e.g. 0, 1, 2, ..., all)
#     - run fgsea exactly as D1_run_fgsea.R does
#       (scoreType = "pos", minSize = 10, maxSize = 500, eps = 0)
#   Then combine species as D1_cross_species_analysis.R does
#     (mean NES across the terms tested in all three species).
#   Record the focal term's rank and percentile, per species and across species.
#
# Nested design over N: within an iteration, each species' neutral background,
# the random order of its focal genes and their selected FST values are drawn
# once. N = k selects the first k genes in that order, so N = 3 contains the
# genes of N = 2, and differences between N reflect only the number of
# selected genes. If a species has fewer focal genes than N, all are selected.
#
# Use s = 0 for a null control: focal genes then get neutral FST too, so the
# focal term's percentile should be uniform (~5% of iterations in the top 5%).
#
# Usage (run from anywhere; paths are absolute or set below):
#   Rscript gsea_power_empirical.R Ne=100k Gen=5 decline=1.00 s=0.1 f0=0.05 \
#       [n_iter=100] [nproc=1] [term=GO:0007596] [out_dir=gsea_empirical_power_nsel]
#       [n_sel=0,1,2,3,4,6,8,10,13,16,20,all]
#
# Cost: each iteration runs 3 species x length(n_sel) fgsea analyses, so the
# default 12 N values cost 12x the old all-genes run. Time a short test first.
#
# Output (one file per condition, appended one iteration at a time, so an
# interrupted job can be resubmitted with the same arguments and continues):
#   <out_dir>/Ne_<Ne>/Gen_<Gen>/d<decline>_s<s>_af<f0>.<term>.nsel.iterations.tsv
#   with a column n_selected (the requested N; capped per species at its
#   number of focal genes, recorded in n_selected_actual)
#
# One row per iteration x species, plus one "cross_species" row per iteration:
#   NES      = mean NES across species (cross_species row)
#   rank_NES = rank by NES (per species) or mean NES (cross_species); 1 = strongest
#   pct_NES  = rank_NES / n_terms * 100 (top x% of terms)
#   rank_p / pct_p = rank by nominal P (per species only)
#   n_species_positive / n_species_sig = species with NES > 0 / padj < 0.05
#                                        (cross_species row only)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

# ============================================================
# Arguments (key=value)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
kv <- strsplit(args, "=", fixed = TRUE)
if (any(lengths(kv) != 2)) stop("Arguments must be key=value, e.g. Ne=100k Gen=5")
opt <- setNames(lapply(kv, `[`, 2), vapply(kv, `[`, "", 1))

get_opt <- function(name, default = NULL) {
  if (!is.null(opt[[name]])) return(opt[[name]])
  if (is.null(default)) stop("Missing required argument: ", name, "=...")
  default
}

NE      <- get_opt("Ne")                       # e.g. 100k (directory spelling)
GEN     <- as.integer(get_opt("Gen"))
DECLINE <- as.numeric(get_opt("decline"))
SEL_S   <- as.numeric(get_opt("s"))            # 0 = null control
F0      <- as.numeric(get_opt("f0"))
N_ITER  <- as.integer(get_opt("n_iter", "100"))
NPROC   <- as.integer(get_opt("nproc", "1"))
TERM    <- get_opt("term", "GO:0007596")       # blood coagulation
OUT_DIR <- get_opt("out_dir", "gsea_empirical_power_nsel")
N_SEL_ARG <- get_opt("n_sel", "0,1,2,3,4,6,8,10,13,16,20,all")

# ============================================================
# Paths
# ============================================================

SIM_DIR  <- get_opt("sim_dir",
  "/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis")
GSEA_DIR <- get_opt("gsea_dir",
  "/xdisk/mcnew/finches/dannyjackson/finches/analyses/gene_set_enrichment")

SPECIES <- c(CRA = "cra", FOR = "for", PAR = "par")

gene_score_files <- file.path(GSEA_DIR, "gene_scores",
                              paste0(SPECIES, ".fst.gene_scores.tsv"))
names(gene_score_files) <- names(SPECIES)
gene_set_file <- file.path(GSEA_DIR, "gene_sets", "go_terms.tsv")

# fgsea settings: identical to D1_run_fgsea.R
FGSEA_MIN_SIZE   <- 10
FGSEA_MAX_SIZE   <- 500
FGSEA_SCORE_TYPE <- "pos"

# ============================================================
# Find simulation files (numeric match on decline / s / f0, as in
# D1_RankPower.top5.test.R, so "1.00" vs "1" spellings don't matter)
# ============================================================

cond_dir <- file.path(SIM_DIR, paste0("Ne_", NE), paste0("Gen_", GEN))
if (!dir.exists(cond_dir)) stop("No simulation directory: ", cond_dir)

sim_files <- list.files(cond_dir, pattern = "\\.summary_stats\\.tsv$", full.names = TRUE)
parsed <- data.table(
  file    = sim_files,
  decline = as.numeric(sub("^d([0-9.]+)_s.*$", "\\1", basename(sim_files))),
  s       = as.numeric(sub("^d[0-9.]+_s([0-9.]+)_af.*$", "\\1", basename(sim_files))),
  af      = as.numeric(sub("^.*_af([0-9.]+)\\.summary_stats\\.tsv$", "\\1", basename(sim_files)))
)

read_fst <- function(file) {
  d <- fread(file, select = "fst_hudson")
  d$fst_hudson[is.finite(d$fst_hudson)]
}

neutral_file <- parsed[abs(decline - DECLINE) < 1e-9 & s == 0 & abs(af - 0.5) < 1e-9]$file[1]
if (is.na(neutral_file)) stop("No neutral file (s = 0, af = 0.5, decline = ", DECLINE, ") in ", cond_dir)
neutral_fst <- read_fst(neutral_file)

if (SEL_S == 0) {
  selected_file <- neutral_file          # null control
  selected_fst  <- neutral_fst
} else {
  selected_file <- parsed[abs(decline - DECLINE) < 1e-9 &
                          abs(s - SEL_S) < 1e-9 & abs(af - F0) < 1e-9]$file[1]
  if (is.na(selected_file)) stop("No selected file for decline = ", DECLINE,
                                 ", s = ", SEL_S, ", f0 = ", F0, " in ", cond_dir)
  selected_fst <- read_fst(selected_file)
}

cat("Neutral file:  ", neutral_file,  " (n =", length(neutral_fst),  ")\n")
cat("Selected file: ", selected_file, " (n =", length(selected_fst), ")\n")
cat(sprintf("Negative FST: neutral %.1f%%, selected %.1f%%\n",
            100 * mean(neutral_fst < 0), 100 * mean(selected_fst < 0)))

# ============================================================
# Empirical gene lists and GO sets (same filtering as D1_run_fgsea.R)
# ============================================================

gene_sets <- unique(fread(gene_set_file)[!is.na(term_id) & !is.na(gene_id),
                                         .(term_id, gene_id)])
pathways  <- split(gene_sets$gene_id, gene_sets$term_id)

if (!TERM %in% names(pathways)) stop(TERM, " not found in ", gene_set_file)

species_genes <- lapply(gene_score_files, function(f) {
  d <- fread(f)[!is.na(mean_score) & !is.na(gene_id)]
  cat(sprintf("%s: %d genes, empirical mean_score range %.4f to %.4f\n",
              basename(f), uniqueN(d$gene_id), min(d$mean_score), max(d$mean_score)))
  unique(d$gene_id)
})

focal_genes <- lapply(species_genes, function(g) intersect(g, pathways[[TERM]]))
for (sp in names(focal_genes)) {
  cat(sprintf("%s: %d of %d %s genes in the ranking\n",
              sp, length(focal_genes[[sp]]), length(pathways[[TERM]]), TERM))
}

if (any(lengths(focal_genes) > length(selected_fst))) {
  stop("More focal genes than selected windows; cannot sample without replacement")
}

# N values: "all" = largest number of focal genes in any species
n_focal_max <- max(lengths(focal_genes))
n_sel_tokens <- trimws(strsplit(N_SEL_ARG, ",", fixed = TRUE)[[1]])
N_SEL <- sort(unique(vapply(n_sel_tokens, function(x)
  if (tolower(x) == "all") n_focal_max else as.integer(x), integer(1))))
N_SEL <- N_SEL[N_SEL >= 0 & N_SEL <= n_focal_max]
if (SEL_S == 0) N_SEL <- 0L    # null control: selected genes would be neutral anyway
cat("N values (selected genes in the term):", paste(N_SEL, collapse = ", "), "\n")

# ============================================================
# Output file and resume
# ============================================================

fmt <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
out_sub <- file.path(OUT_DIR, paste0("Ne_", NE), paste0("Gen_", GEN))
dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_sub, sprintf("d%s_s%s_af%s.%s.nsel.iterations.tsv",
                                       fmt(DECLINE), fmt(SEL_S), fmt(F0),
                                       gsub(":", "", TERM)))

# An (iteration, N) pair counts as done only if its cross-species row was written;
# an iteration is skipped only when all requested N values are done.
done_pairs <- data.table(iteration = integer(0), n_selected = integer(0))
if (file.exists(out_file) && file.info(out_file)$size > 0) {
  prev <- fread(out_file)
  done_pairs <- unique(prev[species == "cross_species", .(iteration, n_selected)])
}
done_iters <- done_pairs[, .(n = uniqueN(n_selected[n_selected %in% N_SEL])), by = iteration][
  n == length(N_SEL), iteration]
cat("Output:", out_file, "| iterations already complete:", length(done_iters), "\n")

# Seed per iteration, unique to the condition, so results are reproducible
# and a resumed job draws the same values it would have drawn originally.
cond_key  <- sprintf("%s_%d_%s_%s_%s_%s", NE, GEN, fmt(DECLINE), fmt(SEL_S), fmt(F0), TERM)
seed_base <- sum(utf8ToInt(cond_key) * seq_along(utf8ToInt(cond_key))) %% 100000 * 1000

# ============================================================
# One fgsea run on simulated scores for one species
# ============================================================

# Draw once per iteration and species: neutral background, a random order of
# the focal genes, and a selected FST for each focal gene (in that order).
draw_species <- function(sp) {
  genes <- species_genes[[sp]]
  # neutral pool (~10,000) is smaller than the gene list, so sample with replacement
  background <- setNames(sample(neutral_fst, length(genes), replace = TRUE), genes)
  fg <- focal_genes[[sp]]
  order_fg <- fg[sample.int(length(fg))]
  sel_vals <- sample(selected_fst, length(fg), replace = FALSE)
  list(background = background, order = order_fg, sel = sel_vals)
}

# Scores with the first n focal genes (in the drawn order) under selection
scores_for_n <- function(draw, n) {
  k <- min(n, length(draw$order))
  scores <- draw$background
  if (k > 0) scores[draw$order[seq_len(k)]] <- draw$sel[seq_len(k)]
  sort(scores, decreasing = TRUE)
}

run_fgsea <- function(stats) {
  res <- fgsea(
    pathways  = pathways,
    stats     = stats,
    minSize   = FGSEA_MIN_SIZE,
    maxSize   = FGSEA_MAX_SIZE,
    scoreType = FGSEA_SCORE_TYPE,
    eps       = 0,
    nproc     = NPROC
  )
  as.data.table(res)[, .(pathway, pval, padj, NES, size)]
}

# Rank and percentile of the focal term within one result table.
# rank 1 = strongest; pct = rank / n_terms * 100 (smaller = stronger)
focal_row <- function(res, iteration, species) {
  res <- res[!is.na(NES)]
  n <- nrow(res)
  res[, rank_NES := frank(-NES, ties.method = "min")]
  res[, rank_p   := frank(pval, ties.method = "min")]
  f <- res[pathway == TERM]
  if (nrow(f) == 0) {
    return(data.table(iteration = iteration, species = species, n_terms = n,
                      NES = NA_real_, pval = NA_real_, padj = NA_real_, size = NA_integer_,
                      rank_NES = NA_integer_, pct_NES = NA_real_,
                      rank_p = NA_integer_, pct_p = NA_real_,
                      n_species_positive = NA_integer_, n_species_sig = NA_integer_))
  }
  data.table(iteration = iteration, species = species, n_terms = n,
             NES = f$NES, pval = f$pval, padj = f$padj, size = f$size,
             rank_NES = f$rank_NES, pct_NES = 100 * f$rank_NES / n,
             rank_p = f$rank_p, pct_p = 100 * f$rank_p / n,
             n_species_positive = NA_integer_, n_species_sig = NA_integer_)
}

# ============================================================
# Iterate
# ============================================================

for (it in setdiff(seq_len(N_ITER), done_iters)) {
  set.seed(seed_base + it)
  draws <- lapply(names(SPECIES), draw_species)
  names(draws) <- names(SPECIES)

  todo_n <- setdiff(N_SEL, done_pairs[iteration == it, n_selected])

 for (n_sel in todo_n) {
  t0 <- Sys.time()
  set.seed(seed_base + it * 1000 + n_sel)   # fgsea permutations; reproducible per (it, N)

  per_species <- lapply(names(SPECIES), function(sp) {
    res <- run_fgsea(scores_for_n(draws[[sp]], n_sel))
    res[, species := sp]
    res
  })
  names(per_species) <- names(SPECIES)

  rows <- rbindlist(lapply(names(SPECIES), function(sp)
    focal_row(per_species[[sp]], it, sp)))

  # ---- Cross-species: as D1_cross_species_analysis.R ----
  wide <- dcast(rbindlist(per_species)[, .(pathway, species, NES)],
                pathway ~ species, value.var = "NES")
  sp_cols <- names(SPECIES)
  wide[, n_species  := rowSums(!is.na(.SD)), .SDcols = sp_cols]
  wide[, n_positive := rowSums(.SD > 0, na.rm = TRUE), .SDcols = sp_cols]
  wide[, mean_NES   := rowMeans(.SD, na.rm = TRUE), .SDcols = sp_cols]
  wide3 <- wide[n_species == 3]
  wide3[, rank_NES := frank(-mean_NES, ties.method = "min")]
  n3 <- nrow(wide3)
  f3 <- wide3[pathway == TERM]

  focal_species <- rows[species %in% sp_cols]
  cross <- data.table(
    iteration = it, species = "cross_species", n_terms = n3,
    NES  = if (nrow(f3)) f3$mean_NES else NA_real_,   # mean NES across species
    pval = NA_real_, padj = NA_real_, size = NA_integer_,
    rank_NES = if (nrow(f3)) f3$rank_NES else NA_integer_,
    pct_NES  = if (nrow(f3)) 100 * f3$rank_NES / n3 else NA_real_,
    rank_p = NA_integer_, pct_p = NA_real_,
    n_species_positive = sum(focal_species$NES > 0, na.rm = TRUE),
    n_species_sig      = sum(focal_species$padj < 0.05, na.rm = TRUE)
  )

  out <- rbind(rows, cross)
  out[, n_selected_actual := c(pmin(n_sel, lengths(focal_genes)[names(SPECIES)]), NA_integer_)]
  out[, `:=`(Ne = NE, Gen = GEN, decline = DECLINE, sel_s = SEL_S, f0 = F0, term = TERM,
             n_selected = n_sel, n_focal_max = n_focal_max)]
  setcolorder(out, c("Ne", "Gen", "decline", "sel_s", "f0", "term", "n_selected"))

  fwrite(out, out_file, sep = "\t", append = file.exists(out_file) &&
           file.info(out_file)$size > 0)

  cat(sprintf("iter %d N=%d: cross-species rank %s of %d (top %.2f%%); per-species pct_NES %s | %.0f s\n",
              it, n_sel, cross$rank_NES, n3, cross$pct_NES,
              paste(sprintf("%s=%.1f", rows$species, rows$pct_NES), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
 }
}

cat("Done:", out_file, "\n")
