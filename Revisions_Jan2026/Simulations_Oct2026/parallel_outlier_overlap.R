#!/usr/bin/env Rscript
# ============================================================
# Parallel-selection overlap test (three species)
#
# Builds on outlier_power_curves.R. Each "species" is a genome of N_GENES
# pseudo-genes (one simulated 10-kb window each) with the SAME gene names in
# all three species, so selection can be assigned to the same genes.
#
# Per selected condition (Ne x Gen x decline x s x f0) and permutation:
#   - every species has N_SEL_PER_SPECIES selected genes, of which
#     N_PARALLEL are selected in all three species (same gene names) and the
#     rest are selected in that species only (different names per species)
#   - neutral genes get FST from the matching neutral file, selected genes
#     from the selected file; each species gets its own independent draws
#   - in each species, take the TOP_N genes by FST
#   - compare the three top lists:
#       overlap3 = genes in the top list of all three species
#       overlap2 = genes in the top list of at least two species
#   - record how many overlap genes are
#       parallel  (selected in all three species)
#       selected  (selected in at least one species)
#       neutral   (selected in none)
#
# With no power, a gene lands in a top-100 list by chance with probability
# TOP_N / N_GENES = 0.01, so the chance 3-way overlap is tiny:
#   expected chance overlap3 = N_GENES * 0.01^3 = 0.01 genes
#   expected chance overlap2 = N_GENES * (3 * 0.01^2 * 0.99 + 0.01^3) ~ 3 genes
#
# Run from the project directory:
#   cd /xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
#   Rscript ${SCRIPTDIR}/parallel_outlier_overlap.R
#
# Writes:
#   output/parallel_overlap.summary.tsv
#   output/parallel_overlap.permutations.tsv
#   output/power_curves_parallel/par<NN>/<metric>/Gen_<g>/
#       power_curves.main.png, power_curves.all_conditions.png,
#       Ne_<NNN>k.decline_<d>.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
})

set.seed(123)

# ============================================================
# Settings
# ============================================================

root_dir <- "."

N_GENES           <- 10000   # pseudo-genes per species (same names in all three)
N_SEL_PER_SPECIES <- 100     # selected genes in each species
PARALLEL_VALUES   <- c(30)   # selected in all three; e.g. c(0, 30, 100) to compare
TOP_N             <- 100     # top list size per species (top 1%)
N_SPECIES         <- 3
N_PERMUTATIONS    <- 100

N_NEUTRAL_PER_SPECIES <- N_GENES - N_SEL_PER_SPECIES

out_dir  <- "output"
out_root <- file.path(out_dir, "power_curves_parallel")

base_Ne      <- "100k"
base_decline <- 1

size_main  <- 13
size_all   <- 10
size_combo <- 12

s_ramp <- c("#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b")

stopifnot(all(PARALLEL_VALUES <= N_SEL_PER_SPECIES),
          N_SPECIES * N_SEL_PER_SPECIES <= N_GENES)

# ============================================================
# Gene names and selection status
# ============================================================

gene_names <- sprintf("gene_%05d", seq_len(N_GENES))

# For a given number of parallel genes, which genes are selected in which species.
# Parallel genes: gene_00001 .. gene_<n_par>, selected in every species.
# Species-specific genes: consecutive non-overlapping blocks after that.
selection_design <- function(n_par) {
  n_specific <- N_SEL_PER_SPECIES - n_par
  sel <- matrix(FALSE, nrow = N_GENES, ncol = N_SPECIES,
                dimnames = list(gene_names, paste0("sp", seq_len(N_SPECIES))))
  if (n_par > 0) sel[seq_len(n_par), ] <- TRUE
  if (n_specific > 0) {
    for (k in seq_len(N_SPECIES)) {
      idx <- n_par + (k - 1) * n_specific + seq_len(n_specific)
      sel[idx, k] <- TRUE
    }
  }
  sel
}

# ============================================================
# One permutation: assign FST, take top lists, score the overlap
# ============================================================

top_genes <- function(fst) {
  shuffle <- sample.int(length(fst))                 # random tie-breaking
  shuffle[order(-fst[shuffle])][seq_len(TOP_N)]
}

one_permutation <- function(null_fst, sel_fst, sel, n_par) {

  # Selected FST: draw all species' selected genes at once without replacement,
  # so no simulated replicate is used twice within a permutation
  sel_draw <- sample(sel_fst, N_SPECIES * N_SEL_PER_SPECIES)

  in_top <- matrix(FALSE, nrow = N_GENES, ncol = N_SPECIES)

  for (k in seq_len(N_SPECIES)) {
    fst <- numeric(N_GENES)
    sel_k <- which(sel[, k])
    neu_k <- which(!sel[, k])
    fst[sel_k] <- sel_draw[(k - 1) * N_SEL_PER_SPECIES + seq_len(N_SEL_PER_SPECIES)]
    fst[neu_k] <- sample(null_fst, N_NEUTRAL_PER_SPECIES)   # without replacement
    in_top[top_genes(fst), k] <- TRUE
  }

  n_top      <- rowSums(in_top)
  n_sel      <- rowSums(sel)
  is_par     <- n_sel == N_SPECIES
  is_any_sel <- n_sel > 0

  o3 <- n_top == N_SPECIES
  o2 <- n_top >= 2

  # Per-species precision for reference (same as outlier_power_curves.R)
  per_species_power <- mean(vapply(seq_len(N_SPECIES),
                                   function(k) mean(sel[in_top[, k], k]),
                                   numeric(1)))

  c(
    overlap3_n         = sum(o3),
    overlap3_parallel  = sum(o3 & is_par),
    overlap3_selected  = sum(o3 & is_any_sel),
    overlap3_neutral   = sum(o3 & !is_any_sel),
    overlap2_n         = sum(o2),
    overlap2_parallel  = sum(o2 & is_par),
    overlap2_selected  = sum(o2 & is_any_sel),
    overlap2_neutral   = sum(o2 & !is_any_sel),
    parallel_in_all_tops = if (n_par > 0) mean(o3[is_par]) else NA_real_,
    per_species_power  = per_species_power
  )
}

read_fst <- function(file) {
  d <- read.delim(file)
  d$fst_hudson[is.finite(d$fst_hudson)]
}

# ============================================================
# Find files (same layout as outlier_power_curves.R)
# ============================================================

all_files <- list.files(root_dir, pattern = "summary_stats\\.tsv$",
                        recursive = TRUE, full.names = TRUE)
all_files <- all_files[str_detect(all_files, "^\\./Ne_[^/]+/Gen_[0-9]+/[^/]+$")]
selected_files <- all_files[!str_detect(basename(all_files), "_s0\\.0_")]

if (length(selected_files) == 0) {
  stop("No selected *.summary_stats.tsv files found under ", root_dir,
       " (run this from the gsea_power_analysis directory)")
}

neutral_cache <- list()

get_neutral <- function(selected_file, decline) {
  cands <- list.files(dirname(selected_file),
                      pattern = "_s0\\.0_af0\\.5\\.summary_stats\\.tsv$",
                      full.names = TRUE)
  if (length(cands) == 0) return(NULL)
  decl <- as.numeric(str_match(basename(cands), "^d([0-9.]+)_s")[, 2])
  hit <- which(decl == as.numeric(decline))
  if (length(hit) == 0) return(NULL)
  f <- cands[hit[1]]
  if (is.null(neutral_cache[[f]])) neutral_cache[[f]] <<- read_fst(f)
  list(file = f, fst = neutral_cache[[f]])
}

# ============================================================
# Run every selected condition
# ============================================================

designs <- setNames(lapply(PARALLEL_VALUES, selection_design), PARALLEL_VALUES)

perm_results <- map_dfr(selected_files, function(selected_file) {

  filename <- basename(selected_file)
  decline  <- str_match(filename, "^d([0-9.]+)_s")[, 2]
  sel_s    <- str_match(filename, "_s([0-9.]+)_af")[, 2]
  af       <- str_match(filename, "_af([0-9.]+)\\.summary")[, 2]
  Ne       <- str_match(selected_file, "Ne_([^/]+)")[, 2]
  Gen      <- str_match(selected_file, "Gen_([0-9]+)")[, 2]

  neutral <- get_neutral(selected_file, decline)
  if (is.null(neutral)) {
    warning("No matching neutral file for: ", selected_file); return(NULL)
  }
  if (length(neutral$fst) < N_NEUTRAL_PER_SPECIES) {
    warning("Not enough usable neutral windows: ", neutral$file,
            " n = ", length(neutral$fst)); return(NULL)
  }
  sel_fst <- read_fst(selected_file)
  if (length(sel_fst) < N_SPECIES * N_SEL_PER_SPECIES) {
    warning("Not enough usable selected windows: ", selected_file,
            " n = ", length(sel_fst)); return(NULL)
  }

  cat("Processing: Ne =", Ne, "Gen =", Gen, "d =", decline,
      "s =", sel_s, "af =", af, "\n")

  map_dfr(PARALLEL_VALUES, function(n_par) {
    sel <- designs[[as.character(n_par)]]
    m <- t(vapply(seq_len(N_PERMUTATIONS),
                  function(p) one_permutation(neutral$fst, sel_fst, sel, n_par),
                  numeric(10)))
    as_tibble(m) %>%
      mutate(Ne = Ne, Gen = as.numeric(Gen), decline = as.numeric(decline),
             sel_s = as.numeric(sel_s), af = as.numeric(af),
             n_parallel = n_par, permutation = seq_len(n()),
             .before = 1)
  })
})

if (nrow(perm_results) == 0) stop("No conditions could be analysed.")

# ============================================================
# Summaries
#   precision3 = parallel genes / all genes in the 3-way overlap, pooled over
#                permutations (NA if the overlap was always empty)
#   recall3    = mean share of the parallel genes found in all three top lists
# ============================================================

pooled <- function(num, den) ifelse(sum(den) > 0, sum(num) / sum(den), NA_real_)

results <- perm_results %>%
  group_by(Ne, Gen, decline, sel_s, af, n_parallel) %>%
  summarise(
    n_permutations         = n(),
    mean_overlap3          = mean(overlap3_n),
    mean_overlap3_parallel = mean(overlap3_parallel),
    mean_overlap3_neutral  = mean(overlap3_neutral),
    prob_overlap3_nonempty = mean(overlap3_n > 0),
    precision3_parallel    = pooled(overlap3_parallel, overlap3_n),
    precision3_selected    = pooled(overlap3_selected, overlap3_n),
    recall3                = mean(parallel_in_all_tops),
    mean_overlap2          = mean(overlap2_n),
    precision2_parallel    = pooled(overlap2_parallel, overlap2_n),
    precision2_selected    = pooled(overlap2_selected, overlap2_n),
    per_species_power      = mean(per_species_power),
    .groups = "drop"
  )

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
write.table(results, file.path(out_dir, "parallel_overlap.summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(perm_results, file.path(out_dir, "parallel_overlap.permutations.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", file.path(out_dir, "parallel_overlap.summary.tsv"))

# ============================================================
# Plots (layout follows outlier_power_curves.R)
#   x = s (log), facets = f0, one figure set per metric x generation time
# ============================================================

ne_num <- function(x) ifelse(grepl("k$", x),
                             as.numeric(sub("k$", "", x)) * 1000, as.numeric(x))

res <- results %>%
  mutate(
    Ne      = as.character(Ne),
    Ne_num  = ne_num(Ne),
    s       = factor(sel_s, levels = sort(unique(sel_s))),
    f0      = factor(paste0("f0 = ", af), levels = paste0("f0 = ", sort(unique(af)))),
    Ne_lab  = factor(paste0("Ne ", Ne), levels = paste0("Ne ", unique(Ne[order(Ne_num)]))),
    Dec_lab = factor(ifelse(decline == 1, "no decline", paste0("decline ", decline)),
                     levels = c("no decline",
                                paste0("decline ",
                                       sort(unique(decline[decline != 1]), decreasing = TRUE))))
  )

s_cols <- setNames(rep_len(s_ramp, nlevels(res$s)), levels(res$s))
s_vals <- sort(unique(res$sel_s))
gens   <- sort(unique(res$Gen))

metrics <- tibble::tribble(
  ~column,               ~tag,           ~label,                                                    ~ymax,
  "precision3_parallel", "precision3",   "Share of 3-way overlap that is parallel selected",        1,
  "recall3",             "recall3",      "Share of parallel genes in all three top lists",          1,
  "mean_overlap3",       "overlap3_n",   "Genes in all three top-100 lists",                        NA,
  "precision2_parallel", "precision2",   "Share of 2+-species overlap that is parallel selected",   1
)

curve_layers <- function(base_size, line_width, ylab, ymax, y_breaks = waiver()) list(
  geom_line(linewidth = line_width, colour = "#3d3d3a"),
  scale_colour_manual(values = s_cols, name = "Selection\ncoefficient (s)", drop = FALSE),
  scale_x_log10(breaks = s_vals, labels = as.character(s_vals)),
  if (is.na(ymax)) scale_y_continuous(limits = c(0, NA), breaks = y_breaks)
  else             scale_y_continuous(limits = c(0, ymax), breaks = y_breaks),
  labs(x = "Selection coefficient (s, log scale)", y = ylab),
  theme_minimal(base_size = base_size),
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
    strip.text       = element_text(face = "bold", colour = "#3d3d3a"),
    strip.text.y     = element_text(angle = 0, hjust = 0),
    axis.text        = element_text(colour = "#6b6a66"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    plot.title       = element_text(face = "bold"),
    legend.position  = "right"
  )
)

for (n_par in PARALLEL_VALUES) {
  for (m in seq_len(nrow(metrics))) {

    col  <- metrics$column[m]
    ylab <- metrics$label[m]
    ymax <- metrics$ymax[m]
    if (n_par == 0 && col %in% c("recall3", "precision3_parallel", "precision2_parallel")) next

    for (g in gens) {
      od <- file.path(out_root, sprintf("par%03d", n_par), metrics$tag[m], paste0("Gen_", g))
      dir.create(od, recursive = TRUE, showWarnings = FALSE)
      dg <- res %>% filter(Gen == g, n_parallel == n_par)
      design_txt <- sprintf("%d selected genes per species, %d selected in all three; top %d of %s genes per species",
                            N_SEL_PER_SPECIES, n_par, TOP_N, format(N_GENES, big.mark = ","))

      d_main <- dg %>% filter(Ne == base_Ne, decline == base_decline)
      if (nrow(d_main) > 0) {
        p <- ggplot(d_main, aes(sel_s, .data[[col]], group = 1)) +
          curve_layers(size_main, 0.8, ylab, ymax) +
          geom_point(aes(colour = s), size = 2) +
          facet_grid(. ~ f0) +
          labs(title = sprintf("Parallel outlier overlap | %s generations", g),
               subtitle = sprintf("Ne = %s, no decline. %s.", base_Ne, design_txt))
        ggsave(file.path(od, "power_curves.main.png"), p,
               width = 14, height = 4.4, dpi = 300, bg = "white")
      }

      n_rows <- n_distinct(dg$Ne_lab) * n_distinct(dg$Dec_lab)
      p_all <- ggplot(dg, aes(sel_s, .data[[col]], group = 1)) +
        curve_layers(size_all, 0.55, ylab, ymax) +
        geom_point(aes(colour = s), size = 1.2) +
        facet_grid(Ne_lab + Dec_lab ~ f0) +
        labs(title = sprintf("Parallel outlier overlap, all conditions | %s generations", g),
             subtitle = design_txt)
      ggsave(file.path(od, "power_curves.all_conditions.png"), p_all,
             width = 14, height = 1.7 * n_rows + 1.6, dpi = 300, bg = "white",
             limitsize = FALSE)

      combos <- dg %>% distinct(Ne, Ne_num, decline, Dec_lab) %>% arrange(Ne_num, desc(decline))
      for (i in seq_len(nrow(combos))) {
        cb <- combos[i, ]
        di <- dg %>% filter(Ne == cb$Ne, decline == cb$decline)
        p_i <- ggplot(di, aes(sel_s, .data[[col]], group = 1)) +
          curve_layers(size_combo, 0.75, ylab, ymax) +
          geom_point(aes(colour = s), size = 1.8) +
          facet_grid(. ~ f0) +
          labs(title = sprintf("Parallel overlap | Ne %s | %s generations | %s",
                               cb$Ne, g, cb$Dec_lab))
        fname <- sprintf("Ne_%sk.decline_%.2f.png",
                         formatC(cb$Ne_num / 1000, width = 3, flag = "0"), cb$decline)
        ggsave(file.path(od, fname), p_i, width = 14, height = 4, dpi = 300, bg = "white")
      }
      message("Saved: ", od)
    }
  }
}

message("\nDone. Plots in ", out_root)