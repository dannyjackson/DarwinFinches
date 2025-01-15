# Tajima's D

# CRA PRE
~/programs/angsd/misc/realSFS ../cra_pre.saf.idx -P 24 > cra_pre.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta ../cra_pre.saf.idx -outname cra_pre -sfs cra_pre.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat cra_pre.thetas.idx

# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../cra_pre.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63039 rows
# 63 sig



max(df_outliers$Tajima)
# 2.679819
min(df_outliers$Tajima)
# -2.012082

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -2.012082
max(outlier_taj_disorder$Tajima)
# -1.783715
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("cra.pre.tajima.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.783715) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.pre.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.pre.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat cra_pre.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 2.530211
min(df_outliers$Tajima)
# -1.96438

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.96438
max(outlier_taj_disorder$Tajima)
# -1.606032
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("cra.pre.tajima.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.606032) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.pre.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.pre.taj.relevantgenenames_50kb.txt


# CRA POST
~/programs/angsd/misc/realSFS ../cra_post.saf.idx -P 24 > cra_post.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta ../cra_post.saf.idx -outname cra_post -sfs cra_post.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat cra_post.thetas.idx


# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat cra_post.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63038 rows
# 63 sig



max(df_outliers$Tajima)
# 3.23515
min(df_outliers$Tajima)
# -1.941603

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.941603
max(outlier_taj_disorder$Tajima)
# -1.532302
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("cra.post.tajima.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.532302) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.post.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' cra.post.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > cra.post.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../cra_post.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 3.19143
min(df_outliers$Tajima)
# -1.601288

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.601288
max(outlier_taj_disorder$Tajima)
# -1.286821
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("cra.post.taj.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.286821) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.post.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' cra.post.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > cra.post.taj.relevantgenenames_50kb.txt






# FOR PRE

~/programs/angsd/misc/realSFS for_pre.saf.idx -P 24 > for_pre.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta for_pre.saf.idx -outname for_pre -sfs for_pre.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat for_pre.thetas.idx


# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../for_pre.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63039 rows
# 63 sig



max(df_outliers$Tajima)
# 3.1093
min(df_outliers$Tajima)
# -1.92363

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.92363
max(outlier_taj_disorder$Tajima)
# -1.354784
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("for.pre.taj.15kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.354784) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.pre.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.pre.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.pre.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../for_pre.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 2.857702
min(df_outliers$Tajima)
# -1.674174

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.674174
max(outlier_taj_disorder$Tajima)
# -1.064564
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("for.pre.taj.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.064564) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.pre.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.pre.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.pre.taj.relevantgenenames_50kb.txt








# FOR POST
~/programs/angsd/misc/realSFS for_pre.saf.idx -P 24 > for_pre.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta for_pre.saf.idx -outname for_pre -sfs for_pre.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat for_pre.thetas.idx


# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../for_post.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63039 rows
# 63 sig



max(df_outliers$Tajima)
# 2.857624
min(df_outliers$Tajima)
# -1.876934

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.876934
max(outlier_taj_disorder$Tajima)
# -1.258204
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("for.post.taj.15kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -1.258204) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.post.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.post.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.post.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../for_post.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 2.657681
min(df_outliers$Tajima)
# -1.52727

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -1.52727
max(outlier_taj_disorder$Tajima)
# -0.968009
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("for.post.taj.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -0.968009) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.post.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.post.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > for.post.taj.relevantgenenames_50kb.txt







# PAR PRE


~/programs/angsd/misc/realSFS par_pre.saf.idx -P 24 > par_pre.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta par_pre.saf.idx -outname par_pre -sfs par_pre.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat par_pre.thetas.idx


# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../par_pre.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63040 rows
# 63 sig



max(df_outliers$Tajima)
# 3.592677
min(df_outliers$Tajima)
# -2.538628

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -2.538628
max(outlier_taj_disorder$Tajima)
# -2.305586
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("par.pre.taj.15kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -2.305586) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.pre.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.pre.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../par_pre.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 3.079961
min(df_outliers$Tajima)
# -2.343898

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -2.343898
max(outlier_taj_disorder$Tajima)
# -2.148685
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("par.pre.taj.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -2.148685) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.pre.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.pre.taj.relevantgenenames_50kb.txt



# PAR POST


~/programs/angsd/misc/realSFS par_post.saf.idx -P 24 > par_post.sfs

# calculate thetas for each site
~/programs/angsd/misc/realSFS saf2theta par_post.saf.idx -outname par_post -sfs par_post.sfs

# estimate tajima's d genome wide
~/programs/angsd/misc/thetaStat do_stat par_post.thetas.idx


# 15 kb windows

# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../par_post.thetas.idx -win 15000 -step 15000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot

sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 63040 rows
# 63 sig



max(df_outliers$Tajima)
# 2.984037
min(df_outliers$Tajima)
# -2.347203

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:63,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -2.347203
max(outlier_taj_disorder$Tajima)
# -2.102294
write_tsv(outlier_taj, "outliertaj_15kb.tsv")




png("par.post.taj.15kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -2.102294) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.post.taj.relevantgenes_15kb.txt


done < outliertaj_15kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.post.taj.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.post.taj.relevantgenenames_15kb.txt



# redo it with 50kb windows
# sliding window tajima's d
~/programs/angsd/misc/thetaStat do_stat ../par_post.thetas.idx -win 50000 -step 50000  -outnames theta.thetasWindow.gz



head -1 theta.thetasWindow.gz.pestPG > theta.thetasWindow.gz.pestPG.chromosomes
grep 'NC_' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> theta.thetasWindow.gz.pestPG.chromosomes
# sed -i 's/NC_//g' theta.thetasWindow.gz.pestPG.chromosomes


echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.foroutliers
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.foroutliers

echo -e 'index\tchr\tmidPos\ttW\ttP\ttF\ttH\ttL\tTajima\tfuf\tfud\tfayh\tzeng\tNsites' > thetasWindow.forplot
grep 'NC' theta.thetasWindow.gz.pestPG | grep -v 'NC_044601.1' >> thetasWindow.forplot
sed -i 's/NC_044571.1/1/g' thetasWindow.forplot
sed -i 's/NC_044572.1/2/g' thetasWindow.forplot
sed -i 's/NC_044573.1/3/g' thetasWindow.forplot
sed -i 's/NC_044574.1/4/g' thetasWindow.forplot
sed -i 's/NC_044575.1/5/g' thetasWindow.forplot
sed -i 's/NC_044576.1/6/g' thetasWindow.forplot
sed -i 's/NC_044577.1/7/g' thetasWindow.forplot
sed -i 's/NC_044578.1/8/g' thetasWindow.forplot
sed -i 's/NC_044579.1/9/g' thetasWindow.forplot
sed -i 's/NC_044580.1/10/g' thetasWindow.forplot
sed -i 's/NC_044581.1/11/g' thetasWindow.forplot
sed -i 's/NC_044582.1/12/g' thetasWindow.forplot
sed -i 's/NC_044583.1/13/g' thetasWindow.forplot
sed -i 's/NC_044584.1/14/g' thetasWindow.forplot
sed -i 's/NC_044585.1/15/g' thetasWindow.forplot
sed -i 's/NC_044586.1/1.1/g' thetasWindow.forplot
sed -i 's/NC_044587.1/17/g' thetasWindow.forplot
sed -i 's/NC_044588.1/18/g' thetasWindow.forplot
sed -i 's/NC_044589.1/19/g' thetasWindow.forplot
sed -i 's/NC_044590.1/20/g' thetasWindow.forplot
sed -i 's/NC_044591.1/21/g' thetasWindow.forplot
sed -i 's/NC_044592.1/22/g' thetasWindow.forplot
sed -i 's/NC_044593.1/23/g' thetasWindow.forplot
sed -i 's/NC_044594.1/24/g' thetasWindow.forplot
sed -i 's/NC_044595.1/25/g' thetasWindow.forplot
sed -i 's/NC_044596.1/26/g' thetasWindow.forplot
sed -i 's/NC_044597.1/27/g' thetasWindow.forplot
sed -i 's/NC_044598.1/28/g' thetasWindow.forplot
sed -i 's/NC_044599.1/29/g' thetasWindow.forplot
sed -i 's/NC_044600.1/4.1/g' thetasWindow.forplot
sed -i 's/NC_044601.1/Z/g' thetasWindow.forplot





# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df_outliers <- read.csv('thetasWindow.foroutliers', sep ='\t')

df <- read.csv('thetasWindow.forplot', sep ='\t')
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




nrow(df_outliers)
# 18882 rows
# 19 sig



max(df_outliers$Tajima)
# 2.567602
min(df_outliers$Tajima)
# -2.314461

ordered_taj <- df_outliers %>% 
 #  orders from smallest to largest
 arrange(Tajima)



outlier_taj_disorder <- ordered_taj[1:19,]

outlier_taj <- outlier_taj_disorder %>% arrange(chr, Tajima)

min(outlier_taj_disorder$Tajima)
# -2.314461
max(outlier_taj_disorder$Tajima)
# -1.987287
write_tsv(outlier_taj, "outliertaj_50kb.tsv")




png("par.post.taj.50kb.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Tajima))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "Tajima's D") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -2.148685) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
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



while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.post.taj.relevantgenes_50kb.txt


done < outliertaj_50kb.tsv 

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' par.post.taj.relevantgenes_50kb.txt | sed 's/ID\=gene\-//g' | sort -u > par.post.taj.relevantgenenames_50kb.txt


# copy gene lists to genelist directory
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/theta/15kb/cra.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/theta/15kb/cra.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/


cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/15kb/for.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/15kb/for.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/


cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/15bk/par.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/15kb/par.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/



# copy plots to copythis
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/theta/15kb/cra.pre.tajima.png /xdisk/mcnew/dannyjackson/copythis/




# copy gene lists to copythis directory
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/theta/15kb/cra.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/theta/15kb/cra.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/


cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/15kb/for.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/15kb/for.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/


cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/15bk/par.pre.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/
cp /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/15kb/par.post.taj.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/