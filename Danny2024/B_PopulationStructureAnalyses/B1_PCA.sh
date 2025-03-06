# B1_PCA.sh with trans, subset
# after evaluating various subsets of the data, I've determined to move forward with this one, which drops one of two genomes that appears to be siblings or close relatives (dropped the one with lower average depth and higher SD of depth)

#!/bin/bash

module load python/3.11/3.11.4
module load R
module load plink/1.9
module load samtools

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all

~/programs/angsd/angsd -vcf-gl genolike.bcf -doPlink 2 -out genolike_plink -doGeno -1 -dopost 1 -domaf 1 -doMajorMinor 1 -sites /xdisk/mcnew/finches/dannyjackson/finches/reference_lists/sites_headless.mafs


plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered

plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 

# bgzip genolike_pruned.beagle.bed

sbatch --account=mcnew \
--job-name=pruneforlinkage \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/pruneforlinkage%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=24:00:00 \
pruneforlinkage.sh

# Submitted batch job 3637791

#!/bin/bash
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all
module load python

# for pop gen
pcangsd -p /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/genolike_pruned -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_likelihoods/all/pruned -t 12 -e 2 --selection --admix


sbatch --account=mcnew \
--job-name=selection_pca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/selection_pca%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=24:00:00 \
selection_pca_subset.sh

Submitted batch job 3643488

# save the following as: /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info_subset.txt

Sample  Species Treatment  
MSB25201    PYRR    Rural
NOCA003 NOCA    Urban
NOCA004 NOCA    Urban
NOCA006 NOCA    Urban
NOCA008 NOCA    Urban
NOCA012 NOCA    Urban
NOCA013 NOCA    Urban
PYRR003 PYRR    Urban
PYRR004 PYRR    Urban
PYRR006 PYRR    Urban
PYRR007 PYRR    Urban
PYRR009 PYRR    Urban
PYRR011 PYRR    Urban
UWBM100619  NOCA    Rural
UWBM100620  NOCA    Rural
UWBM100621  NOCA    Rural
UWBM103345  NOCA    Rural
UWBM103346  PYRR    Rural
UWBM77548   PYRR    Rural
UWBM77718   PYRR    Rural
UWBM77780   PYRR    Rural
UWBM77781   PYRR    Rural
UWBM77856   NOCA    Rural
UWBM77978   NOCA    Rural

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/subset/

C <- as.matrix(read.table("pruned.cov")) # Reads estimated covariance matrix
D <- as.matrix(read.table("output.selection")) # Reads PC based selection statistics
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info_subset.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# Obtain p-values from PC-based selection scan
p <- pchisq(D, 1, lower.tail=FALSE)
write.table(p, "pvalues.txt")

# plot admixture
tbl=read.table("pruned.admix.3.Q")
pdf(file = "admix_pruned.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()





# can we drop one of UWBM77718, UWBM77548?
# drop 718 from this and all future analyses. It has a slightly lower average depth and a higher standard deviation.

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/
echo 'UWBM77718 1' > subsetindv.txt
awk '{print $1}' /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info.txt | tail -n +2 > /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/samplenames.txt

awk 'BEGIN{FS=OFS=" "} NR==FNR{a[NR]=$0; next} {$1=a[FNR]} 1' /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/samplenames.txt genolike_pruned.fam > genolike_pruned.renamed.fam

mv genolike_pruned.renamed.fam genolike_pruned.fam

plink --bed /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.bim --fam /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned.fam --allow-extra-chr --snps-only 'just-acgt' --remove subsetindv.txt --out genolike_pruned_subset --make-bed 


pcangsd -p /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike_pruned_subset -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/pruned_subset -t 12 -e 2 --selection --admix


Sample  Species Treatment  
MSB25201    PYRR    Rural
NOCA003 NOCA    Urban
NOCA004 NOCA    Urban
NOCA006 NOCA    Urban
NOCA008 NOCA    Urban
NOCA012 NOCA    Urban
NOCA013 NOCA    Urban
PYRR003 PYRR    Urban
PYRR004 PYRR    Urban
PYRR006 PYRR    Urban
PYRR007 PYRR    Urban
PYRR009 PYRR    Urban
PYRR011 PYRR    Urban
UWBM100619  NOCA    Rural
UWBM100620  NOCA    Rural
UWBM100621  NOCA    Rural
UWBM103345  NOCA    Rural
UWBM103346  PYRR    Rural
UWBM77548   PYRR    Rural
UWBM77780   PYRR    Rural
UWBM77781   PYRR    Rural
UWBM77856   NOCA    Rural
UWBM77978   NOCA    Rural

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/

C <- as.matrix(read.table("pruned_subset.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists/sample_info_subset.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_subset_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_subset_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_subset_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("pruned_subset.admix.3.Q")
pdf(file = "admix_pruned_subset.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


# plot PCA by each species
# northern cardinals
cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/

plink --bed /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.fam --allow-extra-chr --snps-only 'just-acgt' --remove /xdisk/mcnew/dannyjackson/cardinals/referencelists/subsetnoca.txt --out genolike_pruned_subset_noca --make-bed 

pcangsd -p /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned_subset_noca -o /xdisk/mcnew/dannyjackson/cardinals/analyses/pca/genolike_pruned_subset_noca -t 12 -e 2 --selection --admix

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/pca/

C <- as.matrix(read.table("genolike_pruned_subset_noca.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset_noca.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_noca_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_noca_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_noca_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("genolike_pruned_subset_noca.admix.3.Q")
pdf(file = "admix_pruned_noca.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()

# just pyrrhuloxia

cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/

plink --bed /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned.fam --allow-extra-chr --snps-only 'just-acgt' --remove /xdisk/mcnew/dannyjackson/cardinals/referencelists/subsetpyrr.txt --out genolike_pruned_subset_pyrr --make-bed 

pcangsd -p /xdisk/mcnew/dannyjackson/cardinals/datafiles/geno_likelihoods/all/genolike_pruned_subset_pyrr -o /xdisk/mcnew/dannyjackson/cardinals/analyses/pca/genolike_pruned_subset_pyrr -t 12 -e 2 --selection --admix

# plot it in R
cd /xdisk/mcnew/dannyjackson/cardinals/analyses/pca/

C <- as.matrix(read.table("genolike_pruned_subset_pyrr.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/dannyjackson/cardinals/referencelists/sample_info_subset_pyrr.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$Species)
labs$Treatment <- factor(labs$Treatment)
labs$Sample <- factor(labs$Sample)

pdf(file = "pca_pruned_pyrr_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_pyrr_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_pyrr_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("genolike_pruned_subset_pyrr.admix.3.Q")
pdf(file = "admix_pruned_pyrr.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()







# rerun selection using subset individuals


sbatch --account=mcnew \
--job-name=selection_pca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.selection_pca.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=2:00:00 \
selection_pca_subset.sh

#!/bin/bash

# for selection / this is only able to look at differences between species not urban/rural
module load python

pcangsd -b /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/vcf_likelihoods/all/genolike.beagle.gz -o /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/pca/subset -t 12 -e 2 --selection --pcadapt --sites_save --snp_weights

Submitted batch job 12043363








# done up to here. I did not look further into the PCAdapt output because it is only capable of detecting differences between the species, which is not a major goal of this paper.
