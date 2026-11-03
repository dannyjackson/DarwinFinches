#!/usr/bin/env python3
"""
Same join as Revisions_Jan2026/GSEA/D1_map_GO_finch_chicken.py (finch gene_name == chicken
symbol, inner join), with the file paths as arguments so it can be run on the
propagated table without editing the original script.

Usage:
  python3 map_GO_finch_chicken.py <gene_scores.tsv> <chicken_symbol_GO.tsv> <go_terms_out.tsv>
"""
import sys
import pandas as pd

if len(sys.argv) != 4:
    sys.exit(__doc__)
scores_file, go_file, out_file = sys.argv[1:]

scores = pd.read_csv(scores_file, sep="\t")
go = pd.read_csv(go_file, sep="\t", names=["gene_name", "term_id", "ontology"])
genes = scores[["gene_id", "gene_name"]].drop_duplicates()
merged = genes.merge(go, on="gene_name", how="inner")
out = merged[["term_id", "gene_id", "gene_name"]].drop_duplicates()
out.to_csv(out_file, sep="\t", index=False)

print("Genes in ranking:", genes["gene_id"].nunique())
print("Genes with GO annotation:", merged["gene_id"].nunique())
print("Fraction annotated:", round(merged["gene_id"].nunique() / genes["gene_id"].nunique(), 3))
print("GO terms:", merged["term_id"].nunique())
print("wrote", out_file)
