# Pixy



# DXY but in pixy
#!/bin/bash

#SBATCH --job-name=createvcf_all
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=15:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf_all.%j


cd /xdisk/mcnew/dannyjackson/finches/dxy

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

sbatch createvcf.sh 
Submitted batch job 10129386


module load bcftools 

bcftools convert /xdisk/mcnew/dannyjackson/finches/dxy/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/dxy/allsamples.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/dxy/allsamples.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/dxy/allsamples.vcf

#!/bin/bash

#SBATCH --job-name=createvcf_for
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=15:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf_for.%j


cd /xdisk/mcnew/dannyjackson/finches/dxy/for

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_all_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

sbatch createvcf.sh 
Submitted batch job 2047200


module load bcftools 

bcftools convert angsdput.bcf  -o for.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' for.vcf

sed -i 's/\.realigned\.bam//g' for.vcf

bgzip for.vcf

tabix for.vcf.gz 

#!/bin/bash

#SBATCH --job-name=createvcf_cra
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=15:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf_cra.%j


cd /xdisk/mcnew/dannyjackson/finches/dxy/cra

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_all_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

sbatch createvcf.sh 
Submitted batch job 2047201


module load bcftools 

bcftools convert angsdput.bcf  -o cra.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' cra.vcf

sed -i 's/\.realigned\.bam//g' cra.vcf

bgzip cra.vcf

tabix cra.vcf.gz 

#!/bin/bash

#SBATCH --job-name=createvcf_par
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=15:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf_par.%j


cd /xdisk/mcnew/dannyjackson/finches/dxy/par

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/par_all_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

sbatch createvcf.sh 
Submitted batch job 2047202


module load bcftools 

bcftools convert angsdput.bcf  -o par.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' par.vcf

sed -i 's/\.realigned\.bam//g' par.vcf

bgzip par.vcf

tabix par.vcf.gz 

# FOR pixy

source activate pixy

module load samtools

sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/for_pixy_popfile.txt

#!/bin/bash

#SBATCH --job-name=pixy_for
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_for.%j

pixy --stats pi fst dxy --vcf for.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/for_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/for/pixy --bypass_invariant_check yes

sbatch pixy_for.sh 
Submitted batch job 2047279

# FOR auto pixy

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/for

zcat /xdisk/mcnew/dannyjackson/finches/dxy/for/for.vcf.gz | grep '^#\|^NC_' | grep -v 'NC_044601.1' > for.autosomes.vcf

bgzip for.autosomes.vcf
tabix for.autosomes.vcf.gz

mkdir pixy

#!/bin/bash

#SBATCH --job-name=pixy_for_auto
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=3:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_for_auto.%j


source activate pixy

module load samtools

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/for/pixy

pixy --stats pi fst dxy --vcf ../for.autosomes.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/for_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/for/pixy --bypass_invariant_check yes

sbatch for_auto_pixy.sh 
Submitted batch job 10137072


# CRA pixy
#!/bin/bash

#SBATCH --job-name=pixy_cra
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_cra.%j

source activate pixy

module load samtools

cd /xdisk/mcnew/dannyjackson/finches/dxy/cra/pixy

pixy --stats pi fst dxy --vcf cra.vcf..gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/cra_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/cra/pixy --bypass_invariant_check yes

sbatch pixy_cra.sh 
Submitted batch job 10135229


cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/cra

zcat /xdisk/mcnew/dannyjackson/finches/dxy/cra/cra.vcf.gz | grep '^#\|^NC_' | grep -v 'NC_044601.1' > cra.autosomes.vcf

bgzip cra.autosomes.vcf
tabix cra.autosomes.vcf.gz

mkdir pixy
cd pixy 

#!/bin/bash

#SBATCH --job-name=pixy_cra_auto
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_cra_auto.%j


source activate pixy

module load samtools

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/cra/pixy

pixy --stats pi fst dxy --vcf cra.autosomes.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/cra_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/cra/pixy --bypass_invariant_check yes


sbatch cra_pixy_auto.sh 
Submitted batch job 2047383


# PAR pixy

source activate pixy

module load samtools

sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/par_pixy_popfile.txt

sed -i 's/PARV/lamich_PARV/g' /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/par_pixy_popfile.txt

#!/bin/bash

#SBATCH --job-name=pixy_par
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_par.%j

pixy --stats pi fst dxy --vcf par.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/par_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/par/pixy --bypass_invariant_check yes

sbatch pixy_par.sh 
Submitted batch job 2047286

# pixy par auto

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/par

zcat /xdisk/mcnew/dannyjackson/finches/dxy/par/par.vcf.gz | grep '^#\|^NC_' | grep -v 'NC_044601.1' > par.autosomes.vcf

bgzip par.autosomes.vcf
tabix par.autosomes.vcf.gz

mkdir pixy

#!/bin/bash

#SBATCH --job-name=pixy_par_auto
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_par_auto.%j


source activate pixy

module load samtools

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/par.pixy

pixy --stats pi fst dxy --vcf par.autosomes.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/for_filtering/par_pixy_popfile.txt --window_size 15000 --output_folder /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/par/pixy --bypass_invariant_check yes


sbatch par_auto_pixy.sh 
Submitted batch job 2047384









# Plotting

# CRA autosomes
cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/cra/pixy

head -1 pixy_dxy.txt > cra.auto.pixy_dxy.txt
grep 'NC' pixy_dxy.txt >> cra.auto.pixy_dxy.txt

cp cra.auto.pixy_dxy.txt cra.auto.pixy_dxy.formanhattan.txt

sed -i 's/NC_044571.1/1\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044572.1/2\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044573.1/3\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044574.1/4\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044575.1/5\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044576.1/6\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044577.1/7\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044578.1/8\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044579.1/9\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044580.1/10\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044581.1/11\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044582.1/12\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044583.1/13\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044584.1/14\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044585.1/15\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044586.1/1.1\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044587.1/17\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044588.1/18\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044589.1/19\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044590.1/20\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044591.1/21\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044592.1/22\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044593.1/23\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044594.1/24\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044595.1/25\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044596.1/26\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044597.1/27\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044598.1/28\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044599.1/29\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044600.1/4.1\t/g' cra.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044601.1/999\t/g' cra.auto.pixy_dxy.formanhattan.txt

## dxy 

dxy.all <- read.table("cra.auto.pixy_dxy.txt",header=T)
dxy.subset<-dxy.all[complete.cases(dxy.all),]

SNP<-c(1: (nrow(dxy.subset)))

lower = min(dxy.subset$avg_dxy)
upper = max(dxy.subset$avg_dxy)
cutoff = upper - ((upper-lower)*0.001)

LessThanCutoff <- dxy.subset$avg_dxy < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,dxy.subset)

pdf(file = "cra.auto.dxy_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(dxy.subset$avg_dxy,br=50)
dev.off()

pdf(file = "for.dxy.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(avg_dxy ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], "cra.auto.dxy_sig.csv", row.names=FALSE)



## manhattan

library(qqman)

dxy<-read.table("cra.auto.pixy_dxy.formanhattan.txt", header=TRUE)

dxysubset<-dxy[complete.cases(dxy),]


xu <- mean(dxysubset$avg_dxy)
s <- sd(dxysubset$avg_dxy)
dxysubset$Zdxy = (dxysubset$avg_dxy - xu)/s

SNP<-c(1: (nrow(dxysubset)))

mydf<-data.frame(SNP,dxysubset)


write.csv(mydf[ which(mydf$Zdxy>='5'),], "for.zdxy_sig.csv", row.names=FALSE)


pdf(file = "cra.auto.dxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="avg_dxy",snp="window_pos_1",logp=FALSE,ylab="Dxy"))
dev.off()

pdf(file = "cra.auto.zdxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="Zdxy",snp="window_pos_1",logp=FALSE,ylab="Z DXY", cex = 0.2))
dev.off()




# alternative subest of just the top 0.1%


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('cra.auto.pixy_dxy.txt', sep ='\t')
max(dxy$avg_dxy, na.rm=TRUE)
# 0.6666667
min(dxy$avg_dxy, na.rm=TRUE)
# 0.05555556

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(avg_dxy)) 

nrow(dxy)
# 63097 snps, so top 0.1% would be 63

outlier_dxy_disorder <- ordered_dxy[1:63,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromosome, window_pos_1)

min(outlier_dxy_disorder$avg_dxy)
# 0.5339395
max(outlier_dxy_disorder$avg_dxy)
# 0.6666667
write.csv(outlier_dxy, "outlierdxy.tsv")


# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table("cra.auto.pixy_dxy.formanhattan.txt",header=T)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromosome) %>% 
  summarise(chr_len=max(window_pos_1)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromosome"="chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosome, window_pos_1) %>%
  mutate( BPcum=window_pos_1+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.dxy.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(avg_dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "DXY") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.5339395, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(window_pos_1), alpha=0.7), size=5, force=1.3) +
  
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

cp cra* /xdisk/mcnew/dannyjackson/copythis




cat outlierdxy.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_dxy.tsv

sed -i 's/ /\t/g' outlier_dxy.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $4}' <<<"${line}"`
minpos=$((winstart - 7500))
maxpos=$((winstart + 22500))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.dxy.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_dxy.tsv)




# FOR autosomes

cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/for/pixy

head -1 pixy_dxy.txt > for.auto.pixy_dxy.txt
grep 'NC' pixy_dxy.txt >> for.auto.pixy_dxy.txt

cp for.auto.pixy_dxy.txt for.auto.pixy_dxy.formanhattan.txt

sed -i 's/NC_044571.1/1\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044572.1/2\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044573.1/3\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044574.1/4\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044575.1/5\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044576.1/6\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044577.1/7\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044578.1/8\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044579.1/9\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044580.1/10\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044581.1/11\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044582.1/12\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044583.1/13\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044584.1/14\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044585.1/15\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044586.1/1.1\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044587.1/17\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044588.1/18\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044589.1/19\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044590.1/20\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044591.1/21\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044592.1/22\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044593.1/23\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044594.1/24\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044595.1/25\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044596.1/26\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044597.1/27\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044598.1/28\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044599.1/29\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044600.1/4.1\t/g' for.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044601.1/999\t/g' for.auto.pixy_dxy.formanhattan.txt

## dxy 

dxy.all <- read.table("for.auto.pixy_dxy.txt",header=T)
dxy.subset<-dxy.all[complete.cases(dxy.all),]

SNP<-c(1: (nrow(dxy.subset)))

lower = min(dxy.subset$avg_dxy)
upper = max(dxy.subset$avg_dxy)
cutoff = upper - ((upper-lower)*0.001)

LessThanCutoff <- dxy.subset$avg_dxy < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,dxy.subset)

pdf(file = "for.auto.dxy_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(dxy.subset$avg_dxy,br=50)
dev.off()

pdf(file = "for.dxy.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(avg_dxy ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], "for.auto.dxy_sig.csv", row.names=FALSE)



## manhattan

library(qqman)

dxy<-read.table("for.auto.pixy_dxy.formanhattan.txt", header=TRUE)

dxysubset<-dxy[complete.cases(dxy),]


xu <- mean(dxysubset$avg_dxy)
s <- sd(dxysubset$avg_dxy)
dxysubset$Zdxy = (dxysubset$avg_dxy - xu)/s

SNP<-c(1: (nrow(dxysubset)))

mydf<-data.frame(SNP,dxysubset)


write.csv(mydf[ which(mydf$Zdxy>='5'),], "for.zdxy_sig.csv", row.names=FALSE)


pdf(file = "for.auto.dxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="avg_dxy",snp="window_pos_1",logp=FALSE,ylab="Dxy"))
dev.off()

pdf(file = "for.auto.zdxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="Zdxy",snp="window_pos_1",logp=FALSE,ylab="Z DXY", cex = 0.2))
dev.off()




# alternative subest of just the top 0.1%


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('for.auto.pixy_dxy.txt', sep ='\t')
max(dxy$avg_dxy, na.rm=TRUE)
# 0.5714286
min(dxy$avg_dxy, na.rm=TRUE)
# 0.0625

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(avg_dxy)) 

nrow(dxy)
# 63101 snps, so top 0.1% would be 63

outlier_dxy_disorder <- ordered_dxy[1:63,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromosome, window_pos_1)

min(outlier_dxy_disorder$avg_dxy)
# 0.4765574
max(outlier_dxy_disorder$avg_dxy)
# 0.5714286
write.csv(outlier_dxy, "outlierdxy.tsv")


# for 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table("for.auto.pixy_dxy.txt",header=T)

blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromosome) %>% 
  summarise(chr_len=max(window_pos_1)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromosome"="chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosome, window_pos_1) %>%
  mutate( BPcum=window_pos_1+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.dxy.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(avg_dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "DXY") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.4765574, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(window_pos_1), alpha=0.7), size=5, force=1.3) +
  
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





cat outlierdxy.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_dxy.tsv

sed -i 's/ /\t/g' outlier_dxy.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $4}' <<<"${line}"`
minpos=$((winstart - 7500))
maxpos=$((winstart + 22500))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.dxy.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_dxy.tsv)


cp for* /xdisk/mcnew/dannyjackson/copythis


# PAR autosomes
cd /xdisk/mcnew/dannyjackson/finches/dxy/autosomes/par/pixy

head -1 pixy_dxy.txt > par.auto.pixy_dxy.txt
grep 'NC' pixy_dxy.txt >> par.auto.pixy_dxy.txt

cp par.auto.pixy_dxy.txt par.auto.pixy_dxy.formanhattan.txt

sed -i 's/NC_044571.1/1\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044572.1/2\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044573.1/3\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044574.1/4\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044575.1/5\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044576.1/6\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044577.1/7\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044578.1/8\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044579.1/9\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044580.1/10\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044581.1/11\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044582.1/12\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044583.1/13\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044584.1/14\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044585.1/15\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044586.1/1.1\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044587.1/17\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044588.1/18\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044589.1/19\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044590.1/20\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044591.1/21\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044592.1/22\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044593.1/23\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044594.1/24\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044595.1/25\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044596.1/26\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044597.1/27\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044598.1/28\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044599.1/29\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044600.1/4.1\t/g' par.auto.pixy_dxy.formanhattan.txt
sed -i 's/NC_044601.1/999\t/g' par.auto.pixy_dxy.formanhattan.txt

## dxy 

dxy.all <- read.table("par.auto.pixy_dxy.txt",header=T)
dxy.subset<-dxy.all[complete.cases(dxy.all),]

SNP<-c(1: (nrow(dxy.subset)))

lower = min(dxy.subset$avg_dxy)
upper = max(dxy.subset$avg_dxy)
cutoff = upper - ((upper-lower)*0.001)

LessThanCutoff <- dxy.subset$avg_dxy < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,dxy.subset)

pdf(file = "par.auto.dxy_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(dxy.subset$avg_dxy,br=50)
dev.off()

pdf(file = "for.dxy.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(avg_dxy ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], "par.auto.dxy_sig.csv", row.names=FALSE)



## manhattan

library(qqman)

dxy<-read.table("par.auto.pixy_dxy.formanhattan.txt", header=TRUE)

dxysubset<-dxy[complete.cases(dxy),]


xu <- mean(dxysubset$avg_dxy)
s <- sd(dxysubset$avg_dxy)
dxysubset$Zdxy = (dxysubset$avg_dxy - xu)/s

SNP<-c(1: (nrow(dxysubset)))

mydf<-data.frame(SNP,dxysubset)


write.csv(mydf[ which(mydf$Zdxy>='5'),], "for.zdxy_sig.csv", row.names=FALSE)


pdf(file = "par.auto.dxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="avg_dxy",snp="window_pos_1",logp=FALSE,ylab="Dxy"))
dev.off()

pdf(file = "par.auto.zdxy.manhattan.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="chromosome",bp="window_pos_1",p="Zdxy",snp="window_pos_1",logp=FALSE,ylab="Z DXY", cex = 0.2))
dev.off()




# alternative subest of just the top 0.1%


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('par.auto.pixy_dxy.txt', sep ='\t')
max(dxy$avg_dxy, na.rm=TRUE)
# 0.5538585
min(dxy$avg_dxy, na.rm=TRUE)
# 0.04166667

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(avg_dxy)) 

nrow(dxy)
# 63101 snps, so top 0.1% would be 63

outlier_dxy_disorder <- ordered_dxy[1:63,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromosome, window_pos_1)

min(outlier_dxy_disorder$avg_dxy)
# 0.4826681
max(outlier_dxy_disorder$avg_dxy)
# 0.5538585
write.csv(outlier_dxy, "outlierdxy.tsv")



# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table("par.auto.pixy_dxy.formanhattan.txt",header=T)

blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromosome) %>% 
  summarise(chr_len=max(window_pos_1)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromosome"="chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromosome, window_pos_1) %>%
  mutate( BPcum=window_pos_1+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.dxy.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(avg_dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "DXY") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.4826681, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(window_pos_1), alpha=0.7), size=5, force=1.3) +
  
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







cat outlierdxy.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_dxy.tsv

sed -i 's/ /\t/g' outlier_dxy.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $3}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $4}' <<<"${line}"`
minpos=$((winstart - 7500))
maxpos=$((winstart + 22500))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.dxy.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_dxy.tsv)

cp par* /xdisk/mcnew/dannyjackson/copythis



awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  par.dxy.auto.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  par.dxy.auto.relevantgenenames_15kb.txt

cp par* /xdisk/mcnew/dannyjackson/copythis

rm /xdisk/mcnew/dannyjackson/copythis/*

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  for.dxy.auto.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  for.dxy.auto.relevantgenenames_15kb.txt

cp par.pre.sweed.auto.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/

cp par.post.sweed.auto.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/


