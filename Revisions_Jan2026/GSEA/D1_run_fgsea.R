#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(fgsea)
    library(ggplot2)
    library(GO.db)
    library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(
        paste(
            "Usage:",
            "Rscript D1_run_fgsea.R",
            "<gene_scores.tsv>",
            "<gene_sets.tsv>",
            "<output_prefix>",
            "[score_type]"
        )
    )
}

score_file <- args[1]
gene_set_file <- args[2]
output_prefix <- args[3]

if (length(args) >= 4) {
    score_type <- args[4]
} else {
    score_type <- "auto"
}


# ------------------------------------------------------------
# Read gene scores
# ------------------------------------------------------------

scores_df <- fread(score_file)

required_score_columns <- c(
    "gene_id",
    "mean_score"
)

missing_columns <- setdiff(
    required_score_columns,
    names(scores_df)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing columns in gene score file:",
            paste(missing_columns, collapse = ", ")
        )
    )
}


scores_df <- scores_df[
    !is.na(mean_score) &
    !is.na(gene_id)
]


# If duplicate gene IDs somehow exist, take their mean.
scores_df <- scores_df[
    ,
    .(
        mean_score = mean(mean_score)
    ),
    by = gene_id
]


# ------------------------------------------------------------
# Make ranked vector
# ------------------------------------------------------------

stats <- scores_df$mean_score
names(stats) <- scores_df$gene_id

stats <- sort(
    stats,
    decreasing = TRUE
)


# ------------------------------------------------------------
# Read gene sets
# ------------------------------------------------------------

gene_sets <- fread(gene_set_file)

required_gene_set_columns <- c(
    "term_id",
    "gene_id"
)

missing_columns <- setdiff(
    required_gene_set_columns,
    names(gene_sets)
)

if (length(missing_columns) > 0) {
    stop(
        paste(
            "Missing columns in gene-set file:",
            paste(missing_columns, collapse = ", ")
        )
    )
}


gene_sets <- unique(
    gene_sets[
        !is.na(term_id) &
        !is.na(gene_id),
        .(
            term_id,
            gene_id
        )
    ]
)


# ------------------------------------------------------------
# Report annotation overlap
# ------------------------------------------------------------

annotated_genes <- unique(gene_sets$gene_id)

n_ranked <- length(stats)

n_annotated <- sum(
    names(stats) %in% annotated_genes
)

cat(
    "Genes in ranking:",
    n_ranked,
    "\n"
)

cat(
    "Ranked genes with pathway annotation:",
    n_annotated,
    "\n"
)

cat(
    "Fraction annotated:",
    round(n_annotated / n_ranked, 3),
    "\n"
)


# ------------------------------------------------------------
# Make pathway list
# ------------------------------------------------------------

pathways <- split(
    gene_sets$gene_id,
    gene_sets$term_id
)


# ------------------------------------------------------------
# Determine score type
# ------------------------------------------------------------

if (score_type == "auto") {

    if (all(stats >= 0)) {
        score_type <- "pos"
    } else {
        score_type <- "std"
    }

}


if (!score_type %in% c("pos", "neg", "std")) {
    stop(
        "score_type must be one of: pos, neg, std, auto"
    )
}


cat(
    "GSEA score type:",
    score_type,
    "\n"
)


# ------------------------------------------------------------
# Run GSEA
# ------------------------------------------------------------

set.seed(1)

fgsea_results <- fgsea(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500,
    scoreType = score_type,
    eps = 0
)


fgsea_results <- as.data.table(
    fgsea_results
)


# ------------------------------------------------------------
# Add GO term descriptions
# ------------------------------------------------------------

go_descriptions <- AnnotationDbi::select(
    GO.db,
    keys = unique(fgsea_results$pathway),
    keytype = "GOID",
    columns = c("TERM", "ONTOLOGY")
)

go_descriptions <- as.data.table(go_descriptions)

setnames(
    go_descriptions,
    c("GOID", "TERM", "ONTOLOGY"),
    c("pathway", "term_name", "ontology")
)

# Some GO IDs can return duplicate records; keep one per GO term
go_descriptions <- unique(
    go_descriptions,
    by = "pathway"
)

fgsea_results <- merge(
    fgsea_results,
    go_descriptions,
    by = "pathway",
    all.x = TRUE
)

# If no description is available, fall back to GO ID
fgsea_results[
    is.na(term_name) | term_name == "",
    term_name := pathway
]


# ------------------------------------------------------------
# Convert leading-edge genes to text
# ------------------------------------------------------------

fgsea_results[
    ,
    leadingEdge := vapply(
        leadingEdge,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
    )
]


# Sort by adjusted P value
setorder(
    fgsea_results,
    padj,
    pval
)


# ------------------------------------------------------------
# Write full table
# ------------------------------------------------------------

output_table <- paste0(
    output_prefix,
    ".gsea.tsv"
)

fwrite(
    fgsea_results,
    output_table,
    sep = "\t"
)


# ------------------------------------------------------------
# Write pathways passing FDR < 0.05
# ------------------------------------------------------------

significant <- fgsea_results[
    !is.na(padj) &
    padj < 0.05
]


fwrite(
    significant,
    paste0(
        output_prefix,
        ".gsea.FDR05.tsv"
    ),
    sep = "\t"
)


# ------------------------------------------------------------
# Plot top pathways
# ------------------------------------------------------------

plot_df <- fgsea_results[
    !is.na(NES) &
    !is.na(padj)
]


if (nrow(plot_df) > 0) {

    plot_df <- head(
        plot_df[order(padj)],
        20
    )


    if ("term_name" %in% names(plot_df)) {

        plot_df[
            is.na(term_name) |
            term_name == "",
            term_name := pathway
        ]

        plot_df[
            ,
            label := term_name
        ]

    } else {

        plot_df[
            ,
            label := pathway
        ]

    }


    plot_df[
        ,
        label := factor(
            label,
            levels = rev(label)
        )
    ]


    p <- ggplot(
        plot_df,
        aes(
            x = NES,
            y = label
        )
    ) +
        geom_point(
            aes(
                size = size
            )
        ) +
        theme_bw() +
        labs(
            x = "Normalized enrichment score",
            y = NULL,
            size = "Gene-set size"
        )


    ggsave(
        paste0(
            output_prefix,
            ".top_pathways.pdf"
        ),
        p,
        width = 8,
        height = 6
    )
}


# ------------------------------------------------------------
# Print summary
# ------------------------------------------------------------

cat(
    "Pathways tested:",
    nrow(fgsea_results),
    "\n"
)

cat(
    "Pathways with FDR < 0.05:",
    nrow(significant),
    "\n"
)

cat(
    "Results written to:",
    output_table,
    "\n"
)


# ------------------------------------------------------------
# Plot pathways with strongest positive enrichment
# ------------------------------------------------------------

plot_nes <- fgsea_results[
    !is.na(NES) &
    !is.na(pval)
]

plot_nes <- plot_nes[
    order(-NES)
]

plot_nes <- head(
    plot_nes,
    20
)

if ("term_name" %in% names(plot_nes)) {

    plot_nes[
        is.na(term_name) | term_name == "",
        term_name := pathway
    ]

    plot_nes[, label := term_name]

} else {

    plot_nes[, label := pathway]

}

plot_nes[
    ,
    label := factor(
        label,
        levels = rev(label)
    )
]

p <- ggplot(
    plot_nes,
    aes(
        x = NES,
        y = label
    )
) +
    geom_point(
        aes(
            size = size,
            shape = pval < 0.05
        )
    ) +
    theme_bw() +
    labs(
        x = "Normalized enrichment score (NES)",
        y = NULL,
        size = "Gene-set size",
        shape = "Nominal P < 0.05"
    )

ggsave(
    paste0(
        output_prefix,
        ".top_NES.pdf"
    ),
    p,
    width = 8,
    height = 7
)

plot_df <- fgsea_results[
    !is.na(NES) &
    !is.na(pval) &
    pval > 0
]

plot_df[
    ,
    logp := -log10(pval)
]

p <- ggplot(
    plot_df,
    aes(
        x = NES,
        y = logp
    )
) +
    geom_point() +
    geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed"
    ) +
    theme_bw() +
    labs(
        x = "Normalized enrichment score (NES)",
        y = expression(-log[10](italic(P)))
    )

ggsave(
    paste0(
        output_prefix,
        ".NES_vs_p.pdf"
    ),
    p,
    width = 7,
    height = 6
)




p_hist <- ggplot(
    fgsea_results,
    aes(x = pval)
) +
    geom_histogram(
        bins = 40
    ) +
    theme_bw() +
    labs(
        x = "Nominal P value",
        y = "Number of GO terms"
    )

ggsave(
    paste0(
        output_prefix,
        ".hist.pdf"
    ),
    p_hist,
    width = 7,
    height = 6
)