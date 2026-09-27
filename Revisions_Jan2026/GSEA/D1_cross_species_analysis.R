##### Cross species analyses #####

library(data.table)
library(ggplot2)

cra <- fread("gsea/cra.fst.gsea.tsv")
fortis <- fread("gsea/for.fst.gsea.tsv")
par <- fread("gsea/par.fst.gsea.tsv")


# ------------------------------------------------------------
# Keep pathway descriptions and NES
# ------------------------------------------------------------

cra <- cra[, .(
    pathway,
    term_name,
    CRA = NES
)]

fortis <- fortis[, .(
    pathway,
    term_name,
    FOR = NES
)]

par <- par[, .(
    pathway,
    term_name,
    PAR = NES
)]


# ------------------------------------------------------------
# Merge species
# ------------------------------------------------------------

combined <- Reduce(
    function(x, y)
        merge(
            x,
            y,
            by = c("pathway", "term_name"),
            all = TRUE
        ),
    list(cra, fortis, par)
)


# ------------------------------------------------------------
# Calculate cross-species summary statistics
# ------------------------------------------------------------

combined[
    ,
    mean_NES := rowMeans(
        .SD,
        na.rm = TRUE
    ),
    .SDcols = c("CRA", "FOR", "PAR")
]

combined[
    ,
    n_positive := rowSums(
        .SD > 0,
        na.rm = TRUE
    ),
    .SDcols = c("CRA", "FOR", "PAR")
]

combined[
    ,
    n_species := rowSums(
        !is.na(.SD)
    ),
    .SDcols = c("CRA", "FOR", "PAR")
]


# Keep pathways represented in all three species
combined3 <- combined[
    n_species == 3
]


combined3 <- combined3[
    order(-mean_NES)
]


# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

fwrite(
    combined3,
    "gsea/fst.cross_species.NES.tsv",
    sep = "\t"
)


head(combined3, 30)


# ============================================================
# Plot 1: Top pathways by mean NES across species
# ============================================================

plot_mean <- copy(combined3)

# Focus on pathways enriched positively in all three species
plot_mean <- plot_mean[
    n_positive == 3
]

plot_mean <- head(
    plot_mean[
        order(-mean_NES)
    ],
    25
)


# readable labels
plot_mean[
    ,
    label := paste0(
        term_name,
        " (",
        pathway,
        ")"
    )
]

plot_mean[
    ,
    label := factor(
        label,
        levels = rev(label)
    )
]


p_mean <- ggplot(
    plot_mean,
    aes(
        x = mean_NES,
        y = label
    )
) +
    geom_point(
        size = 3
    ) +
    theme_bw() +
    labs(
        x = "Mean normalized enrichment score",
        y = NULL,
        title = "Pathways positively enriched in all three species"
    ) +
    theme(
        axis.text.y = element_text(size = 8)
    )


ggsave(
    "gsea/fst.cross_species.mean_NES.pdf",
    p_mean,
    width = 10,
    height = 8
)


# ============================================================
# Plot 2: Heatmap of cross-species NES
# ============================================================

heatmap_df <- copy(combined3)

# strongest pathways that are positive in all species
heatmap_df <- heatmap_df[
    n_positive == 3
]

heatmap_df <- head(
    heatmap_df[
        order(-mean_NES)
    ],
    30
)


heatmap_long <- melt(
    heatmap_df,
    id.vars = c(
        "pathway",
        "term_name",
        "mean_NES"
    ),
    measure.vars = c(
        "CRA",
        "FOR",
        "PAR"
    ),
    variable.name = "species",
    value.name = "NES"
)


heatmap_long[
    ,
    label := paste0(
        term_name,
        " (",
        pathway,
        ")"
    )
]


# Order by mean NES
pathway_order <- heatmap_df[
    order(mean_NES),
    paste0(
        term_name,
        " (",
        pathway,
        ")"
    )
]

heatmap_long[
    ,
    label := factor(
        label,
        levels = pathway_order
    )
]


p_heatmap <- ggplot(
    heatmap_long,
    aes(
        x = species,
        y = label,
        fill = NES
    )
) +
    geom_tile() +
    scale_fill_gradient2(
        midpoint = 0
    ) +
    theme_bw() +
    labs(
        x = NULL,
        y = NULL,
        fill = "NES",
        title = "Cross-species pathway enrichment"
    ) +
    theme(
        axis.text.y = element_text(size = 7),
        panel.grid = element_blank()
    )


ggsave(
    "gsea/fst.cross_species.NES_heatmap.pdf",
    p_heatmap,
    width = 9,
    height = 10
)


# ============================================================
# Plot 3: Pairwise species comparisons
# ============================================================

# CRA vs FOR

p_cra_for <- ggplot(
    combined3,
    aes(
        x = CRA,
        y = FOR
    )
) +
    geom_point(
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    theme_bw() +
    labs(
        x = "CRA NES",
        y = "FOR NES",
        title = "CRA vs FOR pathway enrichment"
    )


ggsave(
    "gsea/fst.CRA_vs_FOR.NES.pdf",
    p_cra_for,
    width = 6,
    height = 6
)


# CRA vs PAR

p_cra_par <- ggplot(
    combined3,
    aes(
        x = CRA,
        y = PAR
    )
) +
    geom_point(
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    theme_bw() +
    labs(
        x = "CRA NES",
        y = "PAR NES",
        title = "CRA vs PAR pathway enrichment"
    )


ggsave(
    "gsea/fst.CRA_vs_PAR.NES.pdf",
    p_cra_par,
    width = 6,
    height = 6
)


# FOR vs PAR

p_for_par <- ggplot(
    combined3,
    aes(
        x = FOR,
        y = PAR
    )
) +
    geom_point(
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed"
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dotted"
    ) +
    theme_bw() +
    labs(
        x = "FOR NES",
        y = "PAR NES",
        title = "FOR vs PAR pathway enrichment"
    )


ggsave(
    "gsea/fst.FOR_vs_PAR.NES.pdf",
    p_for_par,
    width = 6,
    height = 6
)


# ============================================================
# Plot 4: NES distributions by species
# ============================================================

all_nes <- melt(
    combined3,
    id.vars = c(
        "pathway",
        "term_name"
    ),
    measure.vars = c(
        "CRA",
        "FOR",
        "PAR"
    ),
    variable.name = "species",
    value.name = "NES"
)


p_dist <- ggplot(
    all_nes,
    aes(
        x = NES,
        fill = species
    )
) +
    geom_density(
        alpha = 0.3
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed"
    ) +
    theme_bw() +
    labs(
        x = "Normalized enrichment score",
        y = "Density",
        fill = "Species",
        title = "Distribution of pathway enrichment across species"
    )


ggsave(
    "gsea/fst.cross_species.NES_distributions.pdf",
    p_dist,
    width = 7,
    height = 5
)
