## CRA
### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/cradiff.txt"


# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $1, $2, $3}' "$TAJIMA" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"


# try a different approach
awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}
END {
    print "chromo\tposition\tfst\ttajima"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA")
    }
}
'  fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" 

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### combine TD and FST
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima")]), ]

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima

df_clean$composite_score <- df_clean$z_fst + df_clean$z_tajima_inv


# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "cra.composite_score.additive.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Get top 1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "cra.composite_score.additive.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# write full table to tsv
write.table(df_clean, "cra.composite_score.additive.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"

# Plot and save to PDF
pdf("cra.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_fst, y = z_tajima_inv)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_FST", y = "z_Tajima_inv", color = "Top 0.1%") +
  theme_minimal()
dev.off()

# Plot and save to PDF
pdf("cra.additive_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = fst, y = tajima)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "FST", y = "Tajima", color = "Top 0.1%") +
  theme_minimal()
dev.off()

```
### Make manhattan plot
#### Prepare file for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thighest_composite' > cra.composite_score.additive.with_chrnum.tsv

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
}' "cra.composite_score.additive.tsv" | tail -n +2 | awk '{print $1, $2, $3, $10, $11}' | tr ' ' '\t'  >> cra.composite_score.additive.with_chrnum.tsv



source ~/programs/DarwinFinches/param_files/cra_params_fst.sh
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
cutoff <- "0.001"  # Convert to numeric
input <- "cra.composite_score.additive.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- as.numeric(min(data$composite_score[data$highest_composite], na.rm = TRUE))

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
  # geom_hline(yintercept = metric_cutoff, linetype="dotted") +
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
# top 0.1 %
awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$10 == "TRUE" { print $1, $2-25000, $2+25000 }' cra.composite_score.additive.tsv | tail -n +2 > cra.composite_score.additive.0.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="cra.composite_score.additive.0.1perc.bed"
GENEFILE="cra.composite_score.additive.0.1perc.genelist.txt"
GENENAMES="cra.composite_score.additive.0.1perc.genenames.txt"
GENEMAPS="cra.composite_score.additive.0.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

# top 1%
awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next } { print $1, $2-25000, $2+25000 }' cra.composite_score.additive.1perc.tsv | tail -n +2 > cra.composite_score.additive.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="cra.composite_score.additive.1perc.bed"
GENEFILE="cra.composite_score.additive.1perc.genelist.txt"
GENENAMES="cra.composite_score.additive.1perc.genenames.txt"
GENEMAPS="cra.composite_score.additive.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

```

## FOR

### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/50000/forpre_forpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/fordiff.txt"

### add chromosome name to Tajima file
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
if [ -n "$CHROM" ]; then
    echo "Adding chromosome conversion column based on conversion file..."

    awk -v convfile="$CHROM" '
    BEGIN {
        FS = OFS = " "  # or use "\t" if you prefer tab-separated
        while ((getline < convfile) > 0) {
            split($0, a, ",")
            conv[a[1]] = a[2]
        }
    }
    FNR == 1 {
        print $0, "new_chromo"
        next
    }
    {
        new_chrom = ($1 in conv ? conv[$1] : $1)
        print $0, new_chrom
    }' "$TAJIMA" | tr ' ' '\t' > "${TAJIMA}.chrom_added.txt"
fi




# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $4, $2, $3}' "${TAJIMA}.chrom_added.txt" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"

# try a different approach
awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}

END {
    print "chromo\tposition\tfst\ttajima"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA")
    }
}
' fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" 

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### combine TD and FST
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima")]), ]

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima

df_clean$composite_score <- df_clean$z_fst + df_clean$z_tajima_inv


# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "for.composite_score.additive.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# Get top 1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "for.composite_score.additive.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# write full table to tsv
write.table(df_clean, "for.composite_score.additive.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"

# Plot and save to PDF
pdf("for.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_fst, y = z_tajima_inv)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_FST", y = "z_Tajima_inv", color = "Top 0.1%") +
  theme_minimal()
dev.off()

# Plot and save to PDF
pdf("for.additive_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = fst, y = tajima)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "FST", y = "Tajima", color = "Top 0.1%") +
  theme_minimal()
dev.off()

```
### Make manhattan plot
#### Prepare file for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thighest_composite' > for.composite_score.additive.with_chrnum.tsv

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
}' "for.composite_score.additive.tsv" | tail -n +2 | awk '{print $1, $2, $3, $10, $11}' | tr ' ' '\t'  >> for.composite_score.additive.with_chrnum.tsv



source ~/programs/DarwinFinches/param_files/for_params_fst.sh

# for some reason this value is causing issues, replace with rounded number:
sed -i 's/10.0276364171942/10.0276/g' for.composite_score.additive.with_chrnum.tsv
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

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for"
color1 <- "#FF817E"
color2 <- "#75002B"
cutoff <- "0.001"  # Convert to numeric
input <- "for.composite_score.additive.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
data <- na.omit(data)

metric_cutoff <- as.numeric(min(data$composite_score[data$highest_composite], na.rm = TRUE))

data$position <- as.numeric(data$position)

data$composite_score <- as.numeric(data$composite_score)

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

plot_data$composite_score=as.numeric(plot_data$composite_score)

plot_data$composite_score <- round(plot_data$composite_score, 4)

cat("Generating plot...\n")
ggplot(plot_data, aes(x = BPcum, y = composite_score)) +
  geom_point(aes(color = as.factor(chromo)), alpha = 0.8, size = 1) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  # geom_hline(yintercept = metric_cutoff, linetype="dotted") +
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
# top 0.1%

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$10 == "TRUE" { print $1, $2-25000, $2+25000 }' for.composite_score.additive.tsv | tail -n +2 > for.composite_score.additive.0.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="for.composite_score.additive.0.1perc.bed"
GENEFILE="for.composite_score.additive.0.1perc.genelist.txt"
GENENAMES="for.composite_score.additive.0.1perc.genenames.txt"
GENEMAPS="for.composite_score.additive.0.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

# top 1%

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
{ print $1, $2-25000, $2+25000 }' for.composite_score.additive.1perc.tsv | tail -n +2 > for.composite_score.additive.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="for.composite_score.additive.1perc.bed"
GENEFILE="for.composite_score.additive.1perc.genelist.txt"
GENENAMES="for.composite_score.additive.1perc.genenames.txt"
GENEMAPS="for.composite_score.additive.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

```
## PAR

### Prepare files
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par

FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/50000/parpre_parpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/pardiff.txt"



### add chromosome name to Tajima file
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
if [ -n "$CHROM" ]; then
    echo "Adding chromosome conversion column based on conversion file..."

    awk -v convfile="$CHROM" '
    BEGIN {
        FS = OFS = " "  # or use "\t" if you prefer tab-separated
        while ((getline < convfile) > 0) {
            split($0, a, ",")
            conv[a[1]] = a[2]
        }
    }
    FNR == 1 {
        print $0, "new_chromo"
        next
    }
    {
        new_chrom = ($1 in conv ? conv[$1] : $1)
        print $0, new_chrom
    }' "$TAJIMA" | tr ' ' '\t' > "${TAJIMA}.chrom_added.txt"
fi




# Temporary sorted files
FST_SORTED="fst.sorted.tmp"
TAJIMA_SORTED="tajima.sorted.tmp"

# Prepare FST: extract chr and midPos and fst
awk 'NR > 1 {print $2, $3, $5}' "$FST" | sort -k1,1 -k2,2n > "$FST_SORTED"

# Prepare TAJIMA: chromo, position, Tajima
awk 'NR > 1 {print $4, $2, $3}' "${TAJIMA}.chrom_added.txt" | sort -k1,1 -k2,2n > "$TAJIMA_SORTED"

# try a different approach
awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    fst[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    tajima[key] = $3
    keys[key] = 1
    next
}

END {
    print "chromo\tposition\tfst\ttajima"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in fst ? fst[k] : "NA"), (k in tajima ? tajima[k] : "NA")
    }
}
' fst.sorted.tmp tajima.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


# Clean up
rm "$FST_SORTED" "$TAJIMA_SORTED" 

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```








### combine TD and FST
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("fst", "tajima")]), ]

### Select only the numeric columns
stat_matrix <- df_clean[, c("fst", "tajima")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))

        fst tajima
fst    1.00   0.03
tajima 0.03   1.00

# Create composite score
### Normalize each component (Z-scores)
df_clean$z_fst <- scale(df_clean$fst)
df_clean$z_tajima <- scale(df_clean$tajima)

### Invert Tajima’s D so that low values contribute positively to the score
df_clean$z_tajima_inv <- -1 * df_clean$z_tajima

df_clean$composite_score <- df_clean$z_fst + df_clean$z_tajima_inv


# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "par.composite_score.additive.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# Get top 1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.99)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "par.composite_score.additive.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


# write full table to tsv
write.table(df_clean, "par.composite_score.additive.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"

# Plot and save to PDF
pdf("par.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_fst, y = z_tajima_inv)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_FST", y = "z_Tajima_inv", color = "Top 0.1%") +
  theme_minimal()
dev.off()

# Plot and save to PDF
pdf("par.additive_composite_score.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = fst, y = tajima)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "FST", y = "Tajima", color = "Top 0.1%") +
  theme_minimal()
dev.off()

```
### Make manhattan plot
#### Prepare file for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thighest_composite' > par.composite_score.additive.with_chrnum.tsv

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
}' "par.composite_score.additive.tsv" | tail -n +2 | awk '{print $1, $2, $3, $10, $11}' | tr ' ' '\t'  >> par.composite_score.additive.with_chrnum.tsv



source ~/programs/DarwinFinches/param_files/par_params_fst.sh
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

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par"
color1 <- "#A6C965"
color2 <- "#203000"
cutoff <- "0.001"  # Convert to numeric
input <- "par.composite_score.additive.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- as.numeric(min(data$composite_score[data$highest_composite], na.rm = TRUE))

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
  # geom_hline(yintercept = metric_cutoff, linetype="dotted") +
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
# top 0.1%

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
$10 == "TRUE" { print $1, $2-25000, $2+25000 }' par.composite_score.additive.tsv | tail -n +2 > par.composite_score.additive.0.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="par.composite_score.additive.0.1perc.bed"
GENEFILE="par.composite_score.additive.0.1perc.genelist.txt"
GENENAMES="par.composite_score.additive.0.1perc.genenames.txt"
GENEMAPS="par.composite_score.additive.0.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

# top 1%

awk 'BEGIN { FS=OFS="\t" }
NR==1 { print "chromo", "position"; next }
 { print $1, $2-25000, $2+25000 }' par.composite_score.additive.1perc.tsv | tail -n +2 > par.composite_score.additive.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="par.composite_score.additive.1perc.bed"
GENEFILE="par.composite_score.additive.1perc.genelist.txt"
GENENAMES="par.composite_score.additive.1perc.genenames.txt"
GENEMAPS="par.composite_score.additive.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

```


## Overlap across species
```
library(ggvenn)

df_cra<-read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra/cra.composite_score.additive.1perc.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for/for.composite_score.additive.1perc.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par/par.composite_score.additive.1perc.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))
# write(intersect_all, file="intersection.fst.all.csv")
# [1] "BMPER" "NUDT4" "TENM3"

union_all <- Reduce(union, list(df_cra, df_for, df_par))
write(union_all, file="union.all.csv")

pdf(file = "intersection.1perc.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

# test for differences
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.txt 
14458 

# cra for
intersect_cra_for <- Reduce(intersect, list(df_cra, df_for))

N=14458 
k = length(intersect_cra_for) # 32
m = length(df_cra) # 488
n = length(df_for) # 321

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))

        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 4.346e-08
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 2.203638 4.861657
sample estimates:
odds ratio 
  3.321546 

# cra par
intersect_cra_par <- Reduce(intersect, list(df_cra, df_par))

N=14458 
k = length(intersect_cra_par) # 30
m = length(df_cra) # 488
n = length(df_par) # 437

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 0.0002227
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.436907 3.208140
sample estimates:
odds ratio 
  2.182747 


# for par
intersect_for_par <- Reduce(intersect, list(df_for, df_par))

N=14458 
k = length(intersect_for_par) # 17
m = length(df_for) # 321
n = length(df_par) # 437

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 0.02951
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.039944 3.008826
sample estimates:
odds ratio 
  1.826352 

p <- c(4.346e-08, 0.0002227, 0.02951)
p.adjust(p, method="fdr")
# 1.3038e-07 3.3405e-04 2.9510e-02

format(p.adjust(p, method="fdr"), scientific = FALSE)
# "0.00000013038" "0.00033405000" "0.02951000000"
# all are sig diff from 0, therefore overlaps are sig

```
## Composite stat across species
```
CRA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra/cra.composite_score.additive.0.1perc.genenames.txt"
FOR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for/for.composite_score.additive.0.1perc.genenames.txt"
PAR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par/par.composite_score.additive.0.1perc.genenames.txt"

CRA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/cra/cra.composite_score.additive.tsv"
FOR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/for/for.composite_score.additive.tsv"
PAR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/par/par.composite_score.additive.tsv"




# Temporary sorted files
CRA_SORTED="cra.sorted.tmp"
FOR_SORTED="for.sorted.tmp"
PAR_SORTED="par.sorted.tmp"

# Prepare files
awk 'NR > 1 {print $1, $2, $9}' "$CRA" | sort -k1,1 -k2,2n > "$CRA_SORTED"
awk 'NR > 1 {print $1, $2, $9}' "$FOR" | sort -k1,1 -k2,2n > "$FOR_SORTED"
awk 'NR > 1 {print $1, $2, $9}' "$PAR" | sort -k1,1 -k2,2n > "$PAR_SORTED"


# try a different approach
awk '
FILENAME == ARGV[1] {
    key = $1 FS $2
    file1[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[2] {
    key = $1 FS $2
    file2[key] = $3
    keys[key] = 1
    next
}
FILENAME == ARGV[3] {
    key = $1 FS $2
    file3[key] = $3
    keys[key] = 1
    next
}
END {
    print "chromo\tposition\tcra\tfor\tpar"
    PROCINFO["sorted_in"] = "@ind_str_asc"
    for (k in keys) {
        split(k, a, FS)
        print a[1], a[2], (k in file1 ? file1[k] : "NA"), (k in file2 ? file2[k] : "NA"), (k in file3 ? file3[k] : "NA")
    }
}
' cra.sorted.tmp for.sorted.tmp par.sorted.tmp | sort -k1,1 -k2,2n | tr ' ' '\t' > combined_stats.tsv


## Clean up
rm "CRA_SORTED" "$FOR_SORTED" "$PAR_SORTED"

wc -l combined_stats.tsv
grep -v 'NA' combined_stats.tsv | wc -l
```




### combine TD and FST
```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("cra", "for.", "par")]), ]

### Select only the numeric columns
stat_matrix <- df_clean[, c("cra", "for.", "par")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))

        cra   for.    par
cra   1.000 -0.004 -0.008
for. -0.004  1.000  0.007
par  -0.008  0.007  1.000

# Create composite score
### Normalize each component (Z-scores)
df_clean$z_cra <- scale(df_clean$cra)
df_clean$z_for <- scale(df_clean$for.)
df_clean$z_par <- scale(df_clean$par)


df_clean$composite_score <- df_clean$z_cra + df_clean$z_for + df_clean$z_par


# Get top 0.1%
df_clean$highest_composite <- "FALSE"
cutoff <- quantile(df_clean$composite_score, 0.999)
top_windows <- df_clean[df_clean$composite_score >= cutoff, ]
df_clean$highest_composite <- df_clean$composite_score >= cutoff
write.table(top_windows, "all.composite_score.additive.0.1perc.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

# write full table to tsv
write.table(df_clean, "all.composite_score.additive.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$highest_composite] <- "CS top 0.1%"

# Plot and save to PDF
pdf("all.cra.for.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_cra, y = z_for)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_cra", y = "z_for", color = "Top 0.1%") +
  theme_minimal()
dev.off()


# Plot and save to PDF
pdf("all.cra.par.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_cra, y = z_par)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_cra", y = "z_par", color = "Top 0.1%") +
  theme_minimal()
dev.off()

# Plot and save to PDF
pdf("all.for.par.additive_composite_score.scaled.pdf", width = 8, height = 6)
ggplot(df_clean, aes(x = z_for, y = z_par)) +
    geom_point(aes(color = highlight_group), size = 1, alpha = 0.8) +
    scale_color_manual(values = c("None" = "gray80",
                                "CS top 0.1%" = "blue")) +
  labs(title = "Plot of Selection Statistics",
       x = "z_for", y = "z_par", color = "Top 0.1%") +
  theme_minimal()
dev.off()

```
### Make manhattan plot
#### Prepare file for plotting
```
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

echo -e 'chromo\tchrom_std\tposition\tcomposite_score\thighest_composite' > all.composite_score.additive.with_chrnum.tsv

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
}' "all.composite_score.additive.tsv" | tail -n +2 | awk '{print $1, $2, $3, $10, $11}' | tr ' ' '\t'  >> all.composite_score.additive.with_chrnum.tsv


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

outdir <- "/xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic/intersection"
color1 <- "#4EAFAF"
color2 <- "#082B64"
cutoff <- "0.001"  # Convert to numeric
input <- "all.composite_score.additive.with_chrnum.tsv"
metric <- "composite_score"


data <- fread(input, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

metric_cutoff <- min(data$composite_score[data$highest_composite], na.rm = TRUE)

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
$10 == "TRUE" { print $1, $2-25000, $2+25000 }' all.composite_score.additive.tsv | tail -n +2 > all.composite_score.additive.0.1perc.bed

GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BEDFILE="all.composite_score.additive.0.1perc.bed"
GENEFILE="all.composite_score.additive.0.1perc.genelist.txt"
GENENAMES="all.composite_score.additive.0.1perc.genenames.txt"
GENEMAPS="all.composite_score.additive.0.1perc.genecoords.txt"

bedtools intersect -a ${GFF} -b ${BEDFILE} -wa > ${GENEFILE}
grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

grep 'ID\=gene' ${GENEFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print($1, $4, $5, arr[1])}' | sed 's/ID\=gene\-//g' | sort -uk4 > ${GENEMAPS}

```




# compute general stats
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas

R

df <-read.csv("crapre/50000/crapre.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)
df <- read.csv("crapost/50000/crapost.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)

df <-read.csv("forpre/50000/forpre.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)
df <- read.csv("forpost/50000/forpost.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)

df <-read.csv("parpre/50000/parpre.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)
df <- read.csv("parpost/50000/parpost.theta.thetasWindow.pestPG.depthmapfiltered", sep="\t", header=FALSE)

colnames(df) <- c("id", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")
mean(df$Tajima)
mean(df$tW/df$nSites)
mean(df$tP/df$nSites)


cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/composite_statistic

### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("cra/combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df <- read.delim("for/combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)
df <- read.delim("par/combined_stats.tsv", header = TRUE, stringsAsFactors = FALSE)

df_clean <- df[complete.cases(df[, c("fst", "tajima")]), ]

mean(df_clean$fst)
mean(df_clean$tajima)
max(df_clean$fst)
min(df_clean$tajima)
max(df_clean$tajima)