import pandas as pd

scores = pd.read_csv(
    "gene_scores/cra.fst.gene_scores.tsv",
    sep="\t"
)

go = pd.read_csv(
    "gene_sets/chicken_symbol_GO_BP.tsv",
    sep="\t",
    names=["gene_name", "term_id", "ontology"]
)

genes = scores[
    ["gene_id", "gene_name"]
].drop_duplicates()

merged = genes.merge(
    go,
    on="gene_name",
    how="inner"
)

out = merged[
    ["term_id", "gene_id", "gene_name"]
].drop_duplicates()

out.to_csv(
    "gene_sets/go_terms.tsv",
    sep="\t",
    index=False
)

print("Genes in ranking:", genes["gene_id"].nunique())
print("Genes with GO annotation:", merged["gene_id"].nunique())
print(
    "Fraction annotated:",
    round(
        merged["gene_id"].nunique() /
        genes["gene_id"].nunique(),
        3
    )
)
print("GO terms:", merged["term_id"].nunique())