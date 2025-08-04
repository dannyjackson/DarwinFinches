#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
GENE <- args[1]
SPECIES <- args[2]

library(ggplot2)
library(reshape2)
library(dplyr)

tab <- read.table(paste0("/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.", SPECIES, ".txt"), header = TRUE)
labs <- data.frame(tab)
labs$Species <- factor(labs$species)
labs$Treatment <- factor(labs$treatment)
labs$Sample <- factor(labs$sample)

C <- as.matrix(read.table(paste0(GENE, ".cov")))
e <- eigen(C)
PC1.PV = round((e$values[1]/sum(e$values))*100, 2)
PC2.PV = round((e$values[2]/sum(e$values))*100, 2)

pca_df <- data.frame(PC1 = e$vectors[,1],
                     PC2 = e$vectors[,2],
                     Species = labs$Species,
                     Treatment = labs$Treatment,
                     Sample = labs$Sample)

plot_pca <- function(df, color_by, outfile) {
  # Calculate convex hulls per group
  hull_df <- df %>%
    group_by(!!sym(color_by)) %>%
    slice(chull(PC1, PC2))

  p <- ggplot(df, aes(x = PC1, y = PC2, color = !!sym(color_by))) +
    geom_point(shape = 16, size = 5) +
    geom_polygon(data = hull_df, aes(x = PC1, y = PC2, fill = !!sym(color_by)), 
                 alpha = 0.2, color = NA, show.legend = FALSE) +
    labs(
      x = paste0("PC1 (", PC1.PV, "%)"),
      y = paste0("PC2 (", PC2.PV, "%)"),
      title = paste0(GENE, " PCA")
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (color_by == "Treatment") {
    p <- p + scale_color_manual(values = c("pre" = "#999999", "post" = "#000000")) +
      scale_fill_manual(values = c("pre" = "#999999", "post" = "#000000"))
  }

  ggsave(outfile, plot = p, width = 7, height = 7)
}

plot_pca(pca_df, "Species", paste0(GENE, "_poly_species.pdf"))
plot_pca(pca_df, "Treatment", paste0(GENE, "_poly_treatment.pdf"))
plot_pca(pca_df, "Sample", paste0(GENE, "_poly_sample.pdf"))
