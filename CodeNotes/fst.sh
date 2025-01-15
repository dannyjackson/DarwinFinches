# Run this code in all three directories:
# /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra
# /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for
# /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par

~/programs/angsd/misc/realSFS fst stats2 here.fst.idx -win 15000 -step 15000 > slidingwindow15kb
~/programs/angsd/misc/realSFS fst stats2 here.fst.idx  -win 50000 -step 50000 >slidingwindow50kb


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_15kb_fst.txt
grep 'NC_' slidingwindow15kb |  grep -v 'NC_044601.1' >> slidingwindow_15kb_fst.txt

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_50kb_fst.txt
grep 'NC_' slidingwindow50kb | grep -v 'NC_044601.1' >> slidingwindow_50kb_fst.txt



sed -i 's/NC_044571.1/1/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_15kb_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_15kb_fst.txt
# sed -i 's/NC_044601.1/Z/g' slidingwindow_15kb_fst.txt


sed -i 's/NC_044571.1/1/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_50kb_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_50kb_fst.txt
# sed -i 's/NC_044601.1/Z/g' slidingwindow_50kb_fst.txt

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_15kb_foroutlier_fst.txt
grep 'NC_' slidingwindow15kb | grep -v 'NC_044601.1' >> slidingwindow_15kb_foroutlier_fst.txt

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_50kb_foroutlier_fst.txt
grep 'NC_' slidingwindow50kb | grep -v 'NC_044601.1'>> slidingwindow_50kb_foroutlier_fst.txt

# CRA

# outliers
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
fst <- read.csv('slidingwindow_15kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.402532
min(fst$fst)
# -0.035135

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

# 62963 snps, so top 0.1% would be 63

outlier_fst_disorder <- ordered_fst[1:63,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.257306
max(outlier_fst_disorder$fst)
# 0.402532
write_tsv(outlier_fst, "outlierfst_15kb.tsv")

## Now with 50kb windows
fst <- read.csv('slidingwindow_50kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.314116
min(fst$fst)
# -0.030907

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

nrow(fst)
# 18876 snps, so top 0.1% would be 19

outlier_fst_disorder <- ordered_fst[1:20,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.204264
max(outlier_fst_disorder$fst)
# 0.314116
write_tsv(outlier_fst, "outlierfst_50kb.tsv")



# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_15kb_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.fst.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.257306, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("cra.fst.15kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.257306, linetype="dashed") +
   
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


### Redo plots with 50kb

df <- read.csv('slidingwindow_50kb_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.fst.50kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.204264, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("cra.fst.50kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.204264, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# Get list of significant genes
echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_15kb_foroutlier_fst.txt
grep 'NC_' slidingwindow | grep -v 'NC_044601.1' >> slidingwindow_15kb_foroutlier_fst.txt






Consider a line in the output file.
Search the gff for genes from the same chromosome.
Only keep genes that are within 15kb of the midpoint of the significant window.
Append those genes to a text file.
Move to the next line in the output file.


# code to find relevant genes from a gff given a list of outlier snps

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_15kb.txt


done < outlierfst_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > relevantgenenames_15kb.txt



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> relevantgenes_50kb.txt


done < outlierfst_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > relevantgenenames_50kb.txt

# FOR


# outliers
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
fst <- read.csv('slidingwindow_15kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.336182
min(fst$fst)
# -0.021022

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

nrow(fst)
# 62974 snps, so top 0.1% would be 63

outlier_fst_disorder <- ordered_fst[1:63,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.181531
max(outlier_fst_disorder$fst)
# 0.336182
write_tsv(outlier_fst, "outlierfst_15kb.tsv")

## Now with 50kb windows
fst <- read.csv('slidingwindow_50kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.245571
min(fst$fst)
# -0.01298

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

nrow(fst)
# 20285 snps, so top 0.1% would be 19

outlier_fst_disorder <- ordered_fst[1:19,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.155116
max(outlier_fst_disorder$fst)
# 0.245571
write_tsv(outlier_fst, "outlierfst_50kb.tsv")



# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_15kb_fst.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.fst.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.181531, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("for.fst.15kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.181531, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


### Redo plots with 50kb

df <- read.csv('slidingwindow_50kb_fst.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.fst.50kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.155116, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("for.fst.50kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.155116, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# Get list of significant genes



Consider a line in the output file.
Search the gff for genes from the same chromosome.
Only keep genes that are within 15kb of the midpoint of the significant window.
Append those genes to a text file.
Move to the next line in the output file.


# code to find relevant genes from a gff given a list of outlier snps

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.relevantgenes_15kb.txt


done < outlierfst_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.relevantgenenames_15kb.txt



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.relevantgenes_50kb.txt


done < outlierfst_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.relevantgenenames_50kb.txt

# PAR


# outliers
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
fst <- read.csv('slidingwindow_15kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.384728
min(fst$fst)
# -0.022581

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

nrow(fst)
# 62981 snps, so top 0.1% would be 63

outlier_fst_disorder <- ordered_fst[1:63,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.210735
max(outlier_fst_disorder$fst)
# 0.384728
write_tsv(outlier_fst, "outlierfst_15kb.tsv")

## Now with 50kb windows
fst <- read.csv('slidingwindow_50kb_foroutlier_fst.txt', sep ='\t')
max(fst$fst)
# 0.265869
min(fst$fst)
# -0.017179

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

nrow(fst)
# 18879 snps, so top 0.1% would be 19

outlier_fst_disorder <- ordered_fst[1:20,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.166938
max(outlier_fst_disorder$fst)
# 0.265869
write_tsv(outlier_fst, "outlierfst_50kb.tsv")



# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_15kb_fst.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.fst.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.210735, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("par.fst.15kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.210735, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


### Redo plots with 50kb

df <- read.csv('slidingwindow_50kb_fst.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.fst.50kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.166938, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





png("par.fst.50kb.scaledyaxis.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$fst))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.166938, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# Get list of significant genes





Consider a line in the output file.
Search the gff for genes from the same chromosome.
Only keep genes that are within 15kb of the midpoint of the significant window.
Append those genes to a text file.
Move to the next line in the output file.


# code to find relevant genes from a gff given a list of outlier snps

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.relevantgenes_15kb.txt


done < outlierfst_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.relevantgenenames_15kb.txt



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.relevantgenes_50kb.txt


done < outlierfst_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.relevantgenenames_50kb.txt


cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/cra.relevantgenes_fst15kb.txt
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/for.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/for.relevantgenes_fst15kb.txt
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/par.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/par.relevantgenes_fst15kb.txt

cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/copythis/cra.relevantgenes_fst15kb.txt

cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/for.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/copythis/for.relevantgenes_fst15kb.txt
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/par.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/copythis/par.relevantgenes_fst15kb.txt

cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/cra.fst.15kb.sigline.png /xdisk/mcnew/dannyjackson/copythis/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/for.fst.15kb.sigline.png /xdisk/mcnew/dannyjackson/copythis/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/par.fst.15kb.sigline.png /xdisk/mcnew/dannyjackson/copythis/