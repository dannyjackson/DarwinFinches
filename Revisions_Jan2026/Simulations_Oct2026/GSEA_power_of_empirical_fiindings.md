# gsea power plus empirical

cd /xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis/gsea_empirical_genes/

```
library(dplyr)
library(readr)

set.seed(123)

# ============================================================
# Input files
# ============================================================

gene_file <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/gene_set_enrichment/gene_scores/cra.fst.gene_scores.tsv"

go_file <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/gene_set_enrichment/gene_sets/chicken_symbol_GO_BP.tsv"

# Change these to the neutral and selected simulation files
neutral_file <- "/path/to/d1.00_s0.0_af0.5.summary_stats.tsv"
selected_file <- "/path/to/d1.00_s0.05_af0.05.summary_stats.tsv"

outfile <- "fst.simulated_blood_coagulation.tsv"

# ============================================================
# Read empirical gene universe
# ============================================================

genes <- read_tsv(gene_file, show_col_types = FALSE) %>%
  select(gene_id, gene_name)

# ============================================================
# Get genes assigned to blood coagulation
# ============================================================

go <- read_tsv(
  go_file,
  col_names = c("gene_name", "GO_term", "ontology"),
  show_col_types = FALSE
)

blood_genes <- go %>%
  filter(GO_term == "GO:0007596") %>%
  distinct(gene_name) %>%
  pull(gene_name)

# ============================================================
# Read simulated FST distributions
# ============================================================

neutral_scores <- read_tsv(
  neutral_file,
  show_col_types = FALSE
) %>%
  filter(is.finite(fst_hudson)) %>%
  pull(fst_hudson)

selected_scores <- read_tsv(
  selected_file,
  show_col_types = FALSE
) %>%
  filter(is.finite(fst_hudson)) %>%
  pull(fst_hudson)

# ============================================================
# Assign simulated scores
# ============================================================

sim_genes <- genes %>%
  mutate(
    # Start by making every gene neutral
    sim_score = sample(
      neutral_scores,
      size = n(),
      replace = TRUE
    ),

    # Mark blood-coagulation genes
    blood_coagulation = gene_name %in% blood_genes
  )

# Replace scores for blood-coagulation genes with selected scores
n_blood <- sum(sim_genes$blood_coagulation)

sim_genes$sim_score[sim_genes$blood_coagulation] <- sample(
  selected_scores,
  size = n_blood,
  replace = TRUE
)

# ============================================================
# Report matching
# ============================================================

cat("Total genes:", nrow(sim_genes), "\n")
cat("Blood coagulation genes in GO file:", length(blood_genes), "\n")
cat("Blood coagulation genes matched to CRA gene file:", n_blood, "\n\n")

cat("Matched blood coagulation genes:\n")
print(
  sim_genes %>%
    filter(blood_coagulation) %>%
    select(gene_id, gene_name, sim_score)
)

cat("\nGO blood-coagulation genes not present in gene-score file:\n")
print(setdiff(blood_genes, genes$gene_name))

# ============================================================
# Save three-column file
# ============================================================

sim_genes %>%
  select(gene_id, gene_name, sim_score) %>%
  write_tsv(outfile)

cat("\nSaved:", outfile, "\n")
```