#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# Read command-line argument
# ------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop(
        "Usage: Rscript goterm_heatmap.R GO:0002250"
    )
}

go_term <- args[1]

# Safe filename version
go_tag <- gsub(":", "", go_term)


# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

go_file <- "gene_sets/go_terms.tsv"

cra_file <- "gene_scores/cra.fst.gene_scores.tsv"
for_file <- "gene_scores/for.fst.gene_scores.tsv"
par_file <- "gene_scores/par.fst.gene_scores.tsv"


# ------------------------------------------------------------
# Read GO annotation
# ------------------------------------------------------------

go <- fread(go_file)

if (!go_term %in% go$term_id) {
    stop(
        paste(
            "GO term not found in",
            go_file,
            ":",
            go_term
        )
    )
}


# ------------------------------------------------------------
# Get genes in requested GO term
# ------------------------------------------------------------

go_genes <- unique(
    go[
        term_id == go_term,
        .(
            gene_id,
            gene_name
        )
    ]
)

cat(
    "GO term:",
    go_term,
    "\n"
)

cat(
    "Genes annotated to GO term:",
    nrow(go_genes),
    "\n"
)


# ------------------------------------------------------------
# Read gene-level FST scores
# ------------------------------------------------------------

cra <- fread(cra_file)[
    ,
    .(
        gene_id,
        gene_name,
        CRA = mean_score
    )
]

fortis <- fread(for_file)[
    ,
    .(
        gene_id,
        gene_name,
        FOR = mean_score
    )
]

par <- fread(par_file)[
    ,
    .(
        gene_id,
        gene_name,
        PAR = mean_score
    )
]


# ------------------------------------------------------------
# Merge GO genes with each species
# ------------------------------------------------------------

cra_go <- merge(
    go_genes,
    cra,
    by = c("gene_id", "gene_name"),
    all.x = TRUE
)

for_go <- merge(
    go_genes,
    fortis,
    by = c("gene_id", "gene_name"),
    all.x = TRUE
)

par_go <- merge(
    go_genes,
    par,
    by = c("gene_id", "gene_name"),
    all.x = TRUE
)


# ------------------------------------------------------------
# Merge species
# ------------------------------------------------------------

combined <- Reduce(
    function(x, y)
        merge(
            x,
            y,
            by = c("gene_id", "gene_name"),
            all = TRUE
        ),
    list(
        cra_go,
        for_go,
        par_go
    )
)


# ------------------------------------------------------------
# Calculate mean FST across species
# ------------------------------------------------------------

combined[
    ,
    mean_FST := rowMeans(
        .SD,
        na.rm = TRUE
    ),
    .SDcols = c(
        "CRA",
        "FOR",
        "PAR"
    )
]


combined <- combined[
    order(-mean_FST)
]


# ------------------------------------------------------------
# Write gene-level table
# ------------------------------------------------------------

output_table <- paste0(
    "gsea/",
    go_tag,
    ".FST.tsv"
)

fwrite(
    combined,
    output_table,
    sep = "\t"
)


# ------------------------------------------------------------
# Convert to long format
# ------------------------------------------------------------

plot_df <- melt(
    combined,
    id.vars = c(
        "gene_id",
        "gene_name",
        "mean_FST"
    ),
    measure.vars = c(
        "CRA",
        "FOR",
        "PAR"
    ),
    variable.name = "species",
    value.name = "FST"
)


# ------------------------------------------------------------
# Order genes by mean FST
# ------------------------------------------------------------

gene_order <- combined[
    order(mean_FST),
    gene_name
]

plot_df[
    ,
    gene_name := factor(
        gene_name,
        levels = gene_order
    )
]


# ------------------------------------------------------------
# Raw FST heatmap with common 0-max scale
# ------------------------------------------------------------

# ------------------------------------------------------------
# Global maximum FST across all genes and species
# ------------------------------------------------------------

max_fst <- max(
    c(
        cra$CRA,
        fortis$FOR,
        par$PAR
    ),
    na.rm = TRUE
)

cat(
    "Maximum FST across all species:",
    max_fst,
    "\n"
)

p_raw <- ggplot(
    plot_df,
    aes(
        x = species,
        y = gene_name,
        fill = FST
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.2
    ) +
    scale_fill_gradient(
        low = "white",
        high = "red",
        limits = c(0, max_fst),
        oob = scales::squish,
        na.value = "grey90"
    ) +
    theme_bw() +
    labs(
        title = go_term,
        subtitle = paste0(
            "Mean FST across windows overlapping each gene; scale = 0 to ",
            round(max_fst, 3)
        ),
        x = NULL,
        y = NULL,
        fill = "Mean FST"
    ) +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10)
    )

ggsave(
    paste0(
        "gsea/",
        go_tag,
        ".FST_heatmap.pdf"
    ),
    p_raw,
    width = 3,
    height = max(
        6,
        nrow(combined) * 0.22
    )
)


# ------------------------------------------------------------
# Row-scaled heatmap
# ------------------------------------------------------------

plot_df[
    ,
    FST_z := {
        if (
            sum(!is.na(FST)) > 1 &&
            sd(FST, na.rm = TRUE) > 0
        ) {
            as.numeric(scale(FST))
        } else {
            NA_real_
        }
    },
    by = gene_name
]


p_scaled <- ggplot(
    plot_df,
    aes(
        x = species,
        y = gene_name,
        fill = FST_z
    )
) +
    geom_tile(
        color = "white",
        linewidth = 0.2
    ) +
    scale_fill_gradient2(
        midpoint = 0,
        na.value = "grey90"
    ) +
    theme_bw() +
    labs(
        title = go_term,
        subtitle = "FST scaled within each gene across species",
        x = NULL,
        y = NULL,
        fill = "Within-gene\nFST z-score"
    ) +
    theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 10)
    )


ggsave(
    paste0(
        "gsea/",
        go_tag,
        ".FST_heatmap_scaled.pdf"
    ),
    p_scaled,
    width = 3,
    height = max(
        6,
        nrow(combined) * 0.22
    )
)


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat(
    "Genes with CRA FST:",
    sum(!is.na(combined$CRA)),
    "\n"
)

cat(
    "Genes with FOR FST:",
    sum(!is.na(combined$FOR)),
    "\n"
)

cat(
    "Genes with PAR FST:",
    sum(!is.na(combined$PAR)),
    "\n"
)

cat(
    "Wrote:",
    output_table,
    "\n"
)

cat(
    "Wrote:",
    paste0(
        "gsea/",
        go_tag,
        ".FST_heatmap.pdf"
    ),
    "\n"
)

cat(
    "Wrote:",
    paste0(
        "gsea/",
        go_tag,
        ".FST_heatmap_scaled.pdf"
    ),
    "\n"
)