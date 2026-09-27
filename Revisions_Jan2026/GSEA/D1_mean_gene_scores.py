#!/usr/bin/env python3

import argparse
import pandas as pd


parser = argparse.ArgumentParser(
    description="""
    Calculate the mean genomic statistic across all windows
    overlapping each gene.

    Expected bedtools intersect input:

    gene_chr
    gene_start
    gene_end
    gene_id
    gene_name
    strand
    window_chr
    window_start
    window_end
    score
    """
)

parser.add_argument(
    "input",
    help="Output from bedtools intersect -wa -wb"
)

parser.add_argument(
    "output",
    help="Output gene score TSV"
)

args = parser.parse_args()


columns = [
    "chrom",
    "gene_start",
    "gene_end",
    "gene_id",
    "gene_name",
    "strand",
    "window_chrom",
    "window_start",
    "window_end",
    "score"
]


df = pd.read_csv(
    args.input,
    sep="\t",
    names=columns,
    header=None
)


df["score"] = pd.to_numeric(
    df["score"],
    errors="coerce"
)

df = df.dropna(subset=["score"])


gene_scores = (
    df
    .groupby(
        [
            "gene_id",
            "gene_name",
            "chrom",
            "gene_start",
            "gene_end"
        ],
        as_index=False
    )
    .agg(
        mean_score=("score", "mean"),
        n_windows=("score", "size")
    )
)


# Rank high values first
gene_scores = gene_scores.sort_values(
    "mean_score",
    ascending=False
)


gene_scores.to_csv(
    args.output,
    sep="\t",
    index=False
)


print(f"Genes scored: {len(gene_scores)}")
print()
print("Number of overlapping windows per gene:")
print(gene_scores["n_windows"].describe())
print()
print(f"Output written to:")
print(args.output)