# effective population size


# Ne from pi and mutation rate

# ~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/dadi/post/all/cra.saf.idx -P 24 > cra.sfs

# Cra pre
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/cra/pre

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > cra.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.saf.idx -sfs cra.sfs -P 24 -outname cra_out 


~/programs/angsd/misc/thetaStat print cra_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat cra_out.thetas.idx

head -1 cra_out.thetas.idx.pestPG > cra_out.thetas.idx.pestPG.chrom
grep 'NC_' cra_out.thetas.idx.pestPG >> cra_out.thetas.idx.pestPG.chrom

awk '{print $5}' cra_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("cra_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))

df[,c("Chr","Ne")]


chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 97,382.73


# Cra post
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > cra.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx -sfs cra.sfs -P 24 -outname cra_out 


~/programs/angsd/misc/thetaStat print cra_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat cra_out.thetas.idx

head -1 cra_out.thetas.idx.pestPG > cra_out.thetas.idx.pestPG.chrom
grep 'NC_' cra_out.thetas.idx.pestPG >> cra_out.thetas.idx.pestPG.chrom

awk '{print $5}' cra_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("cra_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))


chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 89,817.03


# For pre
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/for/pre

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > for.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx -sfs for.sfs -P 24 -outname for_out 


~/programs/angsd/misc/thetaStat print for_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat for_out.thetas.idx

head -1 for_out.thetas.idx.pestPG > for_out.thetas.idx.pestPG.chrom
grep 'NC_' for_out.thetas.idx.pestPG >> for_out.thetas.idx.pestPG.chrom

awk '{print $5}' for_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("for_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))

chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 204,208.4

# For post
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/for/post


~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > for.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx -sfs for.sfs -P 24 -outname for_out 


~/programs/angsd/misc/thetaStat print for_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat for_out.thetas.idx

head -1 for_out.thetas.idx.pestPG > for_out.thetas.idx.pestPG.chrom
grep 'NC_' for_out.thetas.idx.pestPG >> for_out.thetas.idx.pestPG.chrom

awk '{print $5}' for_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("for_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))


chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 215,321.9


# par pre

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/par/pre


~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > par.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.saf.idx -sfs par.sfs -P 24 -outname par_out 


~/programs/angsd/misc/thetaStat print par_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat par_out.thetas.idx

head -1 par_out.thetas.idx.pestPG > par_out.thetas.idx.pestPG.chrom
grep 'NC_' par_out.thetas.idx.pestPG >> par_out.thetas.idx.pestPG.chrom

awk '{print $5}' par_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("par_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))


chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 168,331.2

# par post
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/par/post


~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > par.sfs 

~/programs/angsd/misc/realSFS saf2theta /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx -sfs par.sfs -P 24 -outname par_out 


~/programs/angsd/misc/thetaStat print par_out.thetas.idx 2>/dev/null |head

~/programs/angsd/misc/thetaStat do_stat par_out.thetas.idx

head -1 par_out.thetas.idx.pestPG > par_out.thetas.idx.pestPG.chrom
grep 'NC_' par_out.thetas.idx.pestPG >> par_out.thetas.idx.pestPG.chrom

awk '{print $5}' par_out.thetas.idx.pestPG.chrom

R 
df <- read.csv("par_out.thetas.idx.pestPG.chrom", sep = '\t')
df$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")
df$chrlen <- as.numeric(df$chrlen)

df$wattersonstheta <- df$tW / df$chrlen
# Ne = π/(4μ)
df$wattersonstheta/(4* (2.04*10^-9))
# Mutation rate: 2.04 x 10-9 per base per year

df$Ne <- df$wattersonstheta/(4* (2.04*10^-9))


chrstats <- df[,c("Chr","Ne")]
write.csv(chrstats, "Ne_chrom.csv")

df_auto <- df[-c(31), ]

sum(df_auto$Ne * df_auto$chrlen)/sum(df_auto$chrlen)

Ne = 176974.4

# analyze change over time
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/

R 

df_cra_post <- read.csv("cra/post/Ne_chrom.csv")
df_cra_pre <- read.csv("cra/pre/Ne_chrom.csv")

df_for_post <- read.csv("for/post/Ne_chrom.csv")
df_for_pre <- read.csv("for/pre/Ne_chrom.csv")

df_par_post <- read.csv("par/post/Ne_chrom.csv")
df_par_pre <- read.csv("par/pre/Ne_chrom.csv")

cra_Ne_dff <- df_cra_post$Ne - df_cra_pre$Ne
for_Ne_dff <- df_for_post$Ne - df_for_pre$Ne
par_Ne_dff <- df_par_post$Ne - df_par_pre$Ne

df_change = data.frame(cra_Ne_dff, for_Ne_dff, par_Ne_dff)

df_change$chrlen <- c("115307910","151975198","113180500","71869398","62004293","35429249","37931169","30698694","25429348","20315087","20750786","20242500","18381638","16509717","14057472","70678068","11493000","12127492","11153188","14799631","7903448","4935710","7046529","7720073","3265682","6820980","5764500","5921243","3096680","19696636","70356807")

df_change$chrlen <- as.numeric(df_change$chrlen)

write.csv(df_change, "Ne_chrom.csv")
library(ggplot2)
library(tidyr) 

# Pivot the data into long format
df_long <- pivot_longer(df_change, cols = c(par_Ne_dff, for_Ne_dff, cra_Ne_dff), names_to = "Taxa", values_to = "Ne_value")

# Plot with ggplot2
png("Ne_chr.png", width=2000, height=500)

ggplot(df_long, aes(x = chrlen, y = Ne_value, color = Taxa)) +
  geom_point() +                        # Adds points
  geom_smooth(method = "lm", se = FALSE) + # Adds trend lines
  labs(x = "Chromosome Length", y = "Ne Value", color = "Taxa") +
  theme_minimal()
dev.off()




# cra pre
Ne = 97,382.73

# cra post
Ne = 89,817.03

# for pre
Ne = 204,208.4

# for post
Ne = 215,321.9

# par pre
Ne = 168,331.2

# par post
Ne = 176,974.4

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/Ne/cra

R

df_post <- read.csv("post/Ne_chrom.csv")
df_pre <- read.csv("pre/Ne_chrom.csv")

df_change = df_pre

df_change$Ne_dff <- df_post$Ne - df_pre$Ne

# do chr with highly negative change in Ne also have genes under selection?

min(df_change$Ne_dff)

# chr 20, NC_044590.1, 14799631 len

df_change[order(df_change$Ne_dff), ]


# are any of the blood related genes ID'd in 2 species found on the super negative Ne chroms?
# no, but they're all on negative Ne chroms, just not the most super negative...

GFI1      NC_044578.1 # homologue on NC_044587.1
KIF6      NC_044573.1
STARD13   NC_044571.1
VWF       NC_044578.1 # NC_044586.1
grep 'GFI1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '{print $1}' | sort -u

grep 'KIF6' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '{print $1}' | sort -u

grep 'STARD13' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '{print $1}' | sort -u

grep 'VWF' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '{print $1}' | sort -u


# okay now try to evaluate chromosomal representation of the genes ID'd in CRA
# copied all cra genes from finches_genes.csv to cra_siggenes.txt
sed 's/$/;/g' cra_siggenes.txt > new_cra_siggenes.txt
sed -i 's/^/\\-/g' new_cra_siggenes.txt

grep 'ID=gene' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '{print $9}' |  awk -F';' '{print $1}' | awk -F'-' '{print $2}'

while read -r gene
do 
  grep $gene /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{print $1}' >> file.txt
  
  grep $gene /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{print $9}' |  awk -F';' '{print $1}' | awk -F'-' '{print $2}' >> file.txt

done < new_cra_siggenes.txt

tr '\n' '\t' < file.txt > new_file.txt


sed -i 's/NC/\nNC/g' new_file.txt

grep 'AGA' cra_siggenes.txt
grep 'AGA' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene'




NC_044571.1/1
NC_044572.1/2
NC_044573.1/3
NC_044574.1/4
NC_044575.1/5
NC_044576.1/6
NC_044577.1/7
NC_044578.1/8
NC_044579.1/9
NC_044580.1/10
NC_044581.1/11
NC_044582.1/12
NC_044583.1/13
NC_044584.1/14
NC_044585.1/15
NC_044586.1/1.1
NC_044587.1/17
NC_044588.1/18
NC_044589.1/19
NC_044590.1/20
NC_044591.1/21
NC_044592.1/22
NC_044593.1/23
NC_044594.1/24
NC_044595.1/25
NC_044596.1/26
NC_044597.1/27
NC_044598.1/28
NC_044599.1/29
NC_044600.1/4.1
NC_044601.1/Z






# can i estimate Ne from Neestimator?



sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias\_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/angsty/analyses/raxml/allsamples.vcf


cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/NeEstimator


# making strata
INDIVIDUALS	STRATA
JP4481  CRA_post
JP5410  CRA_post
JP9655  CRA_post
lamich_PARV1  PAR_pre
lamich_PARV2  PAR_pre
lamich_PL15 CRA_pre
lamich_PL16 CRA_pre
lamich_PL4  CRA_pre
lamich_PL7  CRA_pre
lamich_PL9  CRA_pre
RHC097  PAR_post
RHC507  PAR_post
SM031 PAR_post
SM032 PAR_post
SM040 PAR_post
SM059 PAR_post
SM079 PAR_post
SM1067  CRA_post
SM1083  FOR_post
SM1156  FOR_post
SM1157  CRA_post
SM1200  CRA_post
SM1204  FOR_post
SM1231  CRA_post
SM1237  FOR_post
SM1240  CRA_post
SM1266  CRA_post
SM1270  FOR_post
SM1271  FOR_post
SM1272  FOR_post
SM1273  FOR_post
SRR289  FOR_pre
SRR290  FOR_pre
SRR291  FOR_pre
SRR292  FOR_pre
SRR293  FOR_pre
SRR294  FOR_pre
SRR295  FOR_pre
SRR296  FOR_pre
SRR297  FOR_pre
SRR298  FOR_pre
SRR329  PAR_pre
SRR330  PAR_pre
SRR331  PAR_pre
SRR332  PAR_pre
SRR333  PAR_pre
SRR334  PAR_pre
SRR335  PAR_pre
SRR336  PAR_pre
SRR337  PAR_pre
SRR338  PAR_pre

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	
JP4481_all	JP5410_all	JP9655_all	SM1067_all	SM1157_all	SM1200_all	SM1231	SM1240_all	SM1266_all	lamich_PL15	lamich_PL16	lamich_PL4	lamich_PL7	lamich_PL9	SM1083	SM1156	SM1204_all	SM1237	SM1270	SM1271	SM1272	SM1273	SRR2917289	SRR2917290	SRR2917291	SRR2917292	SRR2917293	SRR2917294	SRR2917295	SRR2917296	SRR2917297	SRR2917298	RHC097_all	RHC507_all	SM031_all	SM032_all	SM040_all	SM059_all	SM079_all	lamich_PARV1	lamich_PARV2	SRR2917329	SRR2917330	SRR2917331	SRR2917332	SRR2917333	SRR2917334	SRR2917335	SRR2917336	SRR2917337	SRR2917338
#!/bin/bash

#SBATCH --job-name=transformvcf
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=320gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.transformvcf.%j


conda deactivate > /dev/null 2>&1
IFS=':' read -ra PATHAR <<< "$PATH"
for i in "${PATHAR[@]}"
    do if [[ $i == *"conda"* ]]
        then echo "removing $i from PATH"
    else NEWPATH=$i:$NEWPATH
    fi
done
export PATH=$NEWPATH
module unload gnu8 && module load gnu8
unset NEWPATH
echo "Successfully removed conda"

module load R/4.4.0
# remotes::install_github("thierrygosselin/radiator")

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/NeEstimator

Rscript transformvcf.r

library(radiator)
mystrata <- read.table("pops.txt", header = T)

vcf_data <- read_vcf(data = "/xdisk/mcnew/dannyjackson/finches/angsty/analyses/raxml/allsamples.vcf",
                     strata = mystrata,
                     filter.monomorphic = FALSE,
                     filter.common.markers = FALSE,
                     vcf.stats = FALSE)

## Turn into tidy (should probably be able to do this and previous step in one
## using only this one calling directly to vcf file but it was causing my R
## session to abort session so I do it in two here...maybe not enough RAM)
tidy_data <- tidy_genomic_data(data = vcf_data,
                               strata = mystrata,
                               filter.monomorphic = FALSE,
                               filter.common.markers = FALSE,
                               vcf.stats = FALSE)

## Create chromosome/loci file necessary for NeEstimator to do LDNe method with
## inter-locus comparisons only (by pretending each locus is on its own
## chromosome within using a two column df with matching identical locus and
## chromosome names) across all populations (including all loci)
neestim_chrom_df <- data.frame(CHROM = unique(tidy_data$LOCUS),
                               LOCUS = unique(tidy_data$LOCUS))
write.table(neestim_chrom_df, file = "chrom.txt",
            quote = F, sep = " ", row.names = F, col.names = F)

## For each pop:

for(pop in unique(mystrata$STRATA)) {

idblist <- data.frame(INDIVIDUALS = mystrata[which(mystrata$STRATA != pop),1])

### get correct individuals,
pop_data <- tidy_genomic_data(data = tidy_data,
                              blacklist.id = idblist,
                              filter.monomorphic = FALSE,
                              filter.common.markers = FALSE)

### filter out sites monomorphic in this population,
pop_no_monomorphs <- tidy_genomic_data(data = pop_data,
                                       filter.monomorphic = TRUE,
                                       filter.common.markers = FALSE)

### write genepop file for LDNe in NeEstimator
genomic_converter(data = pop_no_monomorphs,
                  output = "genepop",
                  filename = paste0(as.character(pop),
                                    "_no_monomorphs"))
}

sbatch transform_vcf.sh 
Submitted batch job 10716389
