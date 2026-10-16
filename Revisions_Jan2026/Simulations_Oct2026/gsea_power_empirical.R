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
#     - give every gene in the focal GO term an FST drawn from the SELECTED
#       simulations (all genes in the term are under selection)
#     - run fgsea exactly as D1_run_fgsea.R does
#       (scoreType = "pos", minSize = 10, maxSize = 500, eps = 0)
#   Then combine species as D1_cross_species_analysis.R does
#     (mean NES across the terms tested in all three species).
#   Record the focal term's rank and percentile, per species and across species.
#
# Use s = 0 for a null control: focal genes then get neutral FST too, so the
# focal term's percentile should be uniform (~5% of iterations in the top 5%).
#
# Usage (run from anywhere; paths are absolute or set below):
#   Rscript gsea_power_empirical.R Ne=100k Gen=5 decline=1.00 s=0.1 f0=0.05 \
#       [n_iter=100] [nproc=1] [term=GO:0007596] [out_dir=gsea_empirical_power]
#
# Output (one file per condition, appended one iteration at a time, so an
# interrupted job can be resubmitted with the same arguments and continues):
#   <out_dir>/Ne_<Ne>/Gen_<Gen>/d<decline>_s<s>_af<f0>.<term>.iterations.tsv
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
OUT_DIR <- get_opt("out_dir", "gsea_empirical_power")

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

# ============================================================
# Output file and resume
# ============================================================

fmt <- function(x) format(x, scientific = FALSE, drop0trailing = TRUE)
out_sub <- file.path(OUT_DIR, paste0("Ne_", NE), paste0("Gen_", GEN))
dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_sub, sprintf("d%s_s%s_af%s.%s.iterations.tsv",
                                       fmt(DECLINE), fmt(SEL_S), fmt(F0),
                                       gsub(":", "", TERM)))

done_iters <- integer(0)
if (file.exists(out_file) && file.info(out_file)$size > 0) {
  prev <- fread(out_file)
  # an iteration counts as done only if its cross-species row was written
  done_iters <- unique(prev[species == "cross_species", iteration])
}
cat("Output:", out_file, "| iterations already done:", length(done_iters), "\n")

# Seed per iteration, unique to the condition, so results are reproducible
# and a resumed job draws the same values it would have drawn originally.
cond_key  <- sprintf("%s_%d_%s_%s_%s_%s", NE, GEN, fmt(DECLINE), fmt(SEL_S), fmt(F0), TERM)
seed_base <- sum(utf8ToInt(cond_key) * seq_along(utf8ToInt(cond_key))) %% 100000 * 1000

# ============================================================
# One fgsea run on simulated scores for one species
# ============================================================

simulate_species <- function(sp) {
  genes <- species_genes[[sp]]
  # neutral pool (~10,000) is smaller than the gene list, so sample with replacement
  scores <- setNames(sample(neutral_fst, length(genes), replace = TRUE), genes)
  # every focal gene under selection; independent draws per species
  fg <- focal_genes[[sp]]
  scores[fg] <- sample(selected_fst, length(fg), replace = FALSE)
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
  t0 <- Sys.time()
  set.seed(seed_base + it)

  per_species <- lapply(names(SPECIES), function(sp) {
    res <- run_fgsea(simulate_species(sp))
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
  out[, `:=`(Ne = NE, Gen = GEN, decline = DECLINE, sel_s = SEL_S, f0 = F0, term = TERM)]
  setcolorder(out, c("Ne", "Gen", "decline", "sel_s", "f0", "term"))

  fwrite(out, out_file, sep = "\t", append = file.exists(out_file) &&
           file.info(out_file)$size > 0)

  cat(sprintf("iter %d: cross-species rank %s of %d (top %.2f%%); per-species pct_NES %s | %.0f s\n",
              it, cross$rank_NES, n3, cross$pct_NES,
              paste(sprintf("%s=%.1f", rows$species, rows$pct_NES), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

cat("Done:", out_file, "\n")
