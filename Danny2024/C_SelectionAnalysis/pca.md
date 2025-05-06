# Compute PCA of sliding window statistics
## CRA
### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/cradiff.txt"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/cra/deltaAF_lrt_full.windowedavg.tsv.header"



# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"
DELTA_AF_SORTED="delta_af.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $1, $2, $3}' "$TAJIMA" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"

# Prepare DELTA_AF: chromo, position, delta_af
awk 'NR > 1 {print $1, $2+25000, $4}' "$DELTAAF" | sort -k1,1 -k2,2n > "$DELTA_AF_SORTED"

# try a different approach
awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    delta[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[3] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}
END {
    print "chromo\tposition\tfst\ttajima\tdelta_af"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA"), (k in delta ? delta[k] : "NA")
    }
}
' delta_af.sorted.tmp fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" "$DELTA_AF_SORTED"

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### compute the actual PCA
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima", "delta_af")]), ]
df_clean$delta_af <- abs(df_clean$delta_af)

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima", "delta_af")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


           fst tajima delta_af
fst      1.000  0.098    0.355
tajima   0.098  1.000    0.144
delta_af 0.355  0.144    1.000


### Plot the three stats in 3D space
library(plotly)
library(htmlwidgets)

# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_daf <- scale(df_clean$delta_af)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima


# Perform PCA on z-scored stats
pca_input <- df_clean[, c("z_fst", "z_daf", "z_tajima_inv")]
pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)


# Append PCA scores to dataframe
df_clean$PC1 <- pca_result$x[, 1]
df_clean$PC2 <- pca_result$x[, 2]

### Loadings = variable contribution directions
loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale

# compute composite score
## first, extract variances explained by PC1 and PC2
pca_var <- summary(pca_result)$importance
pc1_var <- pca_var["Proportion of Variance", "PC1"]
pc2_var <- pca_var["Proportion of Variance", "PC2"]

## quick visual of PC to see which direction of PCs we need to sample

# Plot and save to PDF
pdf("pca_plot.quicklook.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = PC1, y = PC2)) +
    geom_point(color = "gray", size = 1, alpha = 0.8) +
  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +
  labs(title = "PCA of Selection Statistics",
       x = "PC1", y = "PC2", color = "Top 1%") +
  theme_minimal()
dev.off()

### we want positive PC1 and negative PC2
df_clean$composite_score <-  df_clean$PC1 * (-1 * df_clean$PC2)

# Get top 1%
df_clean$high_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$high_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "cra.composite_score.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "cra.composite_score.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# write full table to tsv
write.table(df_clean, "cra.composite_score.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Label top windows based on composite stat
df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$high_composite] <- "CS top 1%"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"

# Plot and save to PDF
pdf("pca_plot_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 1%" = "red",
                                "CS top 0.1%" = "blue")) +
  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +
  labs(title = "PCA of Selection Statistics",
       x = "PC1", y = "PC2", color = "Top 1%") +
  theme_minimal()
dev.off()



# plot in 3D space with z-stats
# Define the color map (matching desired categories)
color_map <- c("None" = "gray", "CS top 1%" = "red", "CS top 0.1%" = "blue")

# Ensure highlight_group is a factor with these levels
df_clean$highlight_group <- factor(df_clean$highlight_group,
                                    levels = names(color_map))

# 3D plot
p <- plot_ly(
  data = df_clean,
  x = ~fst,
  y = ~tajima,
  z = ~delta_af,
  type = "scatter3d",
  mode = "markers",
  color = ~highlight_group,
  colors = color_map,
  marker = list(size = 2, opacity = 0.8)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Z FST"),
      yaxis = list(title = "Z inverse Tajima's D"),
      zaxis = list(title = "Z |ΔAF|")
    ),
    legend = list(title = list(text = "Highlight Group"))
  )

# Save as HTML
saveWidget(p, "cra.composite_score_3Dplot.html", selfcontained = FALSE)



```
### Make manhattan plot
#### Prepare file for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thigh_composite\thighest_composite' > cra.composite_score.with_chrnum.tsv

awk -F'\t' 'BEGIN {
    FS=OFS="\t"
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[2]] = a[1]
    }
}
NR==1 {
    print "chromo", $0
    next
}
{
    print map[$1], $0
}' "cra.composite_score.tsv" | tail -n +2 | awk '{print $1, $2, $3, $13, $14, $15}' | tr ' ' '\t'  >> cra.composite_score.with_chrnum.tsv



source ~/programs/DarwinFinches/paraam_files/cra_params_fst.sh
```
#### Plot manhattan plot
```
# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://forn.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra"
color1 <- "#4EAFAF"
color2 <- "#082B64"
cutoff <- "0.01"  # Convert to numeric
input <- "cra.composite_score.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- min(data$composite_score[data$high_composite], na.rm = TRUE)
metric_cutoff2 <- min(data$composite_score[data$highest_composite], na.rm = TRUE)

data$position <- as.numeric(data$position)

data[[metric]] <- as.numeric(data[[metric]])

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chromo <- factor(data$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")
ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  geom_hline(yintercept = metric_cutoff2, linetype="dotted") +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0(input, ".sigline.png"), 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
```
### Filter windows and create list of genes
```
awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$13 == "TRUE" { print $1, $2-25000, $2+25000 }' cra.composite_score.tsv | tail -n +2 > cra.composite_score.high.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="cra.composite_score.high.bed"
GENEFILE="cra.composite_score.high.genelist.txt"
GENENAMES="cra.composite_score.high.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$14 == "TRUE" { print $1, $2-25000, $2+25000 }' cra.composite_score.tsv | tail -n +2 > cra.composite_score.highest.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="cra.composite_score.highest.bed"
GENEFILE="cra.composite_score.highest.genelist.txt"
GENENAMES="cra.composite_score.highest.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}
```




## FOR

### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/50000/forpre_forpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/fordiff.txt"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/for/deltaAF_lrt_full.windowedavg.tsv.header"

# Add a column to the Tajima file with the 

awk -F'\t' 'BEGIN {
    FS=OFS=" "
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[1]] = a[2]
    }
}
NR==1 {
    print "chrom_std", $0
    next
}
{
    print map[$1], $0
}' "$TAJIMA" | tr ' ' '\t' >> tajima_chroms.tsv

TAJIMA="tajima_chroms.tsv"


# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"
DELTA_AF_SORTED="delta_af.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $1, $3, $4}' "$TAJIMA" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"

# Prepare DELTA_AF: chromo, position, delta_af
awk 'NR > 1 {print $1, $2+25000, $4}' "$DELTAAF" | sort -k1,1 -k2,2n > "$DELTA_AF_SORTED"

awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    delta[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[3] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}
END {
    print "chromo\tposition\tfst\ttajima\tdelta_af"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA"), (k in delta ? delta[k] : "NA")
    }
}
' delta_af.sorted.tmp fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" "$DELTA_AF_SORTED"

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### compute the actual PCA
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima", "delta_af")]), ]
df_clean$delta_af <- abs(df_clean$delta_af)

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima", "delta_af")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


          fst tajima delta_af
fst      1.00  0.060    0.470
tajima   0.06  1.000    0.027
delta_af 0.47  0.027    1.000

### Plot the three stats in 3D space
library(plotly)
library(htmlwidgets)

# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_daf <- scale(df_clean$delta_af)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima

### Combine into a composite score
df_clean$composite_score <- df_clean$z_fst * df_clean$z_daf * df_clean$z_tajima_inv

# Rank-based score
df_clean$rank_fst <- rank(df_clean$fst) / nrow(df_clean)
df_clean$rank_daf <- rank(df_clean$delta_af) / nrow(df_clean)
df_clean$rank_tajima <- 1 - (rank(df_clean$tajima) / nrow(df_clean))  # low D ranks higher

df_clean$rank_score <- df_clean$rank_fst * df_clean$rank_daf * df_clean$rank_tajima

# Get top 1%
df_clean$high_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$high_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "for.top_windows_composite_score.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "for.top_windows_composite_score.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# write full table to tsv
write.table(df_clean, "for.composite_score.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Perform PCA on z-scored stats
pca_input <- df_clean[, c("z_fst", "z_daf", "z_tajima_inv")]
pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)

# Append PCA scores to dataframe
df_clean$PC1 <- pca_result$x[, 1]
df_clean$PC2 <- pca_result$x[, 2]

# Label top windows
df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$high_composite] <- "CS top 1%"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"


### Loadings = variable contribution directions
loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale


# Plot and save to PDF
pdf("pca_plot_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 1%" = "red",
                                "CS top 0.1%" = "blue")) +
  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", size = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +
  labs(title = "PCA of Selection Statistics",
       x = "PC1", y = "PC2", color = "Top 1%") +
  theme_minimal()
dev.off()



# plot in 3D space with z-stats
# Define the color map (matching desired categories)
color_map <- c("None" = "gray", "CS top 1%" = "red", "CS top 0.1%" = "blue")

# Ensure highlight_group is a factor with these levels
df_clean$highlight_group <- factor(df_clean$highlight_group,
                                    levels = names(color_map))

# 3D plot
p <- plot_ly(
  data = df_clean,
  x = ~fst,
  y = ~tajima,
  z = ~delta_af,
  type = "scatter3d",
  mode = "markers",
  color = ~highlight_group,
  colors = color_map,
  marker = list(size = 2, opacity = 0.8)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "FST"),
      yaxis = list(title = "Tajima's D"),
      zaxis = list(title = "|ΔAF|")
    ),
    legend = list(title = list(text = "Highlight Group"))
  )

# Save as HTML
saveWidget(p, "for.composite_score_3Dplot.html", selfcontained = FALSE)



```
### Make manhattan plot
#### Prepare the data for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thigh_composite\thighest_composite' > for.composite_score.with_chrnum.tsv

awk -F'\t' 'BEGIN {
    FS=OFS="\t"
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[2]] = a[1]
    }
}
NR==1 {
    print "chromo", $0
    next
}
{
    print map[$1], $0
}' "for.composite_score.tsv" | tail -n +2 | awk '{print $1, $2, $3, $11, $16, $17}' | tr ' ' '\t'  >> for.composite_score.with_chrnum.tsv



source ~/programs/DarwinFinches/param_files/for_params_fst.sh
```
#### Plot the manhattan plot
```

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://forn.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for"
color1 <- "#FF817E"
color2 <- "#75002B"
cutoff <- "0.01"  # Convert to numeric
input <- "for.composite_score.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- min(data$composite_score[data$high_composite], na.rm = TRUE)
metric_cutoff2 <- min(data$composite_score[data$highest_composite], na.rm = TRUE)

data$position <- as.numeric(data$position)

data[[metric]] <- as.numeric(data[[metric]])

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chromo <- factor(data$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")
ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  geom_hline(yintercept = metric_cutoff2, linetype="dotted") +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0(input, ".sigline.png"), 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
```
### Filter windows and create list of genes
```
awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$15 == "TRUE" { print $1, $2-25000, $2+25000 }' for.composite_score.tsv | tail -n +2 > for.composite_score.high.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="for.composite_score.high.bed"
GENEFILE="for.composite_score.high.genelist.txt"
GENENAMES="for.composite_score.high.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$16 == "TRUE" { print $1, $2-25000, $2+25000 }' for.composite_score.tsv | tail -n +2 > for.composite_score.highest.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="for.composite_score.highest.bed"
GENEFILE="for.composite_score.highest.genelist.txt"
GENENAMES="for.composite_score.highest.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}
```



## PAR

### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/50000/parpre_parpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/pardiff.txt"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/par/deltaAF_lrt_full.windowedavg.tsv.header"


# Add a column to the Tajima file with the 

awk -F'\t' 'BEGIN {
    FS=OFS=" "
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[1]] = a[2]
    }
}
NR==1 {
    print "chrom_std", $0
    next
}
{
    print map[$1], $0
}' "$TAJIMA" | tr ' ' '\t' >> tajima_chroms.tsv

TAJIMA="tajima_chroms.tsv"


# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"
DELTA_AF_SORTED="delta_af.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $1, $3, $4}' "$TAJIMA" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"

# Prepare DELTA_AF: chromo, position, delta_af
awk 'NR > 1 {print $1, $2+25000, $4}' "$DELTAAF" | sort -k1,1 -k2,2n > "$DELTA_AF_SORTED"

awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    delta[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[3] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}
END {
    print "chromo\tposition\tfst\ttajima\tdelta_af"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA"), (k in delta ? delta[k] : "NA")
    }
}
' delta_af.sorted.tmp fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" "$DELTA_AF_SORTED"

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### compute the actual PCA
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima", "delta_af")]), ]
df_clean$delta_af <- abs(df_clean$delta_af)

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima", "delta_af")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


           fst tajima delta_af
fst      1.000  0.031    0.459
tajima   0.031  1.000   -0.004
delta_af 0.459 -0.004    1.000

### Plot the three stats in 3D space
library(plotly)
library(htmlwidgets)

# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_daf <- scale(df_clean$delta_af)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima

### Combine into a composite score
df_clean$composite_score <- df_clean$z_fst * df_clean$z_daf * df_clean$z_tajima_inv

# Rank-based score
df_clean$rank_fst <- rank(df_clean$fst) / nrow(df_clean)
df_clean$rank_daf <- rank(df_clean$delta_af) / nrow(df_clean)
df_clean$rank_tajima <- 1 - (rank(df_clean$tajima) / nrow(df_clean))  # low D ranks higher

df_clean$rank_score <- df_clean$rank_fst * df_clean$rank_daf * df_clean$rank_tajima

# Get top 1%
df_clean$high_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$high_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "par.top_windows_composite_score.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "par.top_windows_composite_score.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# write full table to tsv
write.table(df_clean, "par.composite_score.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Perform PCA on z-scored stats
pca_input <- df_clean[, c("z_fst", "z_daf", "z_tajima_inv")]
pca_result <- prcomp(pca_input, center = TRUE, scale. = FALSE)

# Append PCA scores to dataframe
df_clean$PC1 <- pca_result$x[, 1]
df_clean$PC2 <- pca_result$x[, 2]

# Label top windows
df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$high_composite] <- "CS top 1%"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"


### Loadings = variable contribution directions
loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale


# Plot and save to PDF
pdf("pca_plot_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 1%" = "red",
                                "CS top 0.1%" = "blue")) +
  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", size = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +
  labs(title = "PCA of Selection Statistics",
       x = "PC1", y = "PC2", color = "Top 1%") +
  theme_minimal()
dev.off()



# plot in 3D space with z-stats
# Define the color map (matching desired categories)
color_map <- c("None" = "gray", "CS top 1%" = "red", "CS top 0.1%" = "blue")

# Ensure highlight_group is a factor with these levels
df_clean$highlight_group <- factor(df_clean$highlight_group,
                                    levels = names(color_map))

# 3D plot
p <- plot_ly(
  data = df_clean,
  x = ~fst,
  y = ~tajima,
  z = ~delta_af,
  type = "scatter3d",
  mode = "markers",
  color = ~highlight_group,
  colors = color_map,
  marker = list(size = 2, opacity = 0.8)
) %>%
  layout(
    scene = list(
      xaxis = list(title = "FST"),
      yaxis = list(title = "Tajima's D"),
      zaxis = list(title = "|ΔAF|")
    ),
    legend = list(title = list(text = "Highlight Group"))
  )

# Save as HTML
saveWidget(p, "par.composite_score_3Dplot.html", selfcontained = FALSE)



```
# make manhattan plot
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thigh_composite\thighest_composite' > par.composite_score.with_chrnum.tsv

awk -F'\t' 'BEGIN {
    FS=OFS="\t"
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[2]] = a[1]
    }
}
NR==1 {
    print "chromo", $0
    next
}
{
    print map[$1], $0
}' "par.composite_score.tsv" | tail -n +2 | awk '{print $1, $2, $3, $11, $16, $17}' | tr ' ' '\t'  >> par.composite_score.with_chrnum.tsv



source ~/programs/DarwinFinches/param_files/par_params_fst.sh


# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://forn.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par"
color1 <- "#A6C965"
color2 <- "#203000"
cutoff <- "0.01"  # Convert to numeric
input <- "par.composite_score.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- min(data$composite_score[data$high_composite], na.rm = TRUE)
metric_cutoff2 <- min(data$composite_score[data$highest_composite], na.rm = TRUE)

data$position <- as.numeric(data$position)

data[[metric]] <- as.numeric(data[[metric]])

# Prepare data for plotting
cat("Preparing data for plotting...\n")
data$chromo <- factor(data$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data <- data %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))


# Plot
cat("Generating plot...\n")
ggplot(plot_data, aes(x = BPcum, y = !!sym(metric))) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  geom_hline(yintercept = metric_cutoff2, linetype="dotted") +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

ggsave(filename = paste0(input, ".sigline.png"), 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")
```
### Filter windows and create list of genes
```
awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$15 == "TRUE" { print $1, $2-25000, $2+25000 }' par.composite_score.tsv | tail -n +2 > par.composite_score.high.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="par.composite_score.high.bed"
GENEFILE="par.composite_score.high.genelist.txt"
GENENAMES="par.composite_score.high.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$16 == "TRUE" { print $1, $2-25000, $2+25000 }' par.composite_score.tsv | tail -n +2 > par.composite_score.highest.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="par.composite_score.highest.bed"
GENEFILE="par.composite_score.highest.genelist.txt"
GENENAMES="par.composite_score.highest.genenames.txt"
bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}




































