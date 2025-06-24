# B1_PCA.sh with trans, subset

#!/bin/bash

module load python/3.11/3.11.4
module load R
module load plink/1.9
module load samtools

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt -anc /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -doPlink 2 -out genolike_plink -doGeno -1 -dopost 1 


sbatch --account=mcnew \
--job-name=angsd_pca \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/angsd_pca%j \
--nodes=1 \
--ntasks-per-node=4 \
--time=24:00:00 \
angsd_pca.sh


plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered

plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 

# bgzip genolike_pruned.beagle.bed



# Submitted batch job 3637791

#!/bin/bash
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all
module load python
module load bcftools

# for pop gen
pcangsd -p /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/pruned -t 12 -e 2 --selection --admix

# save the following as: /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt 

sample  species treatment
JP4481  CRA     post
JP5410  CRA     post
JP9655  CRA     post
SM1067  CRA     post
SM1157  CRA     post
SM1200  CRA     post
SM1231  CRA     post
SM1240  CRA     post
SM1266  CRA     post
lamich_PL15     CRA     pre
lamich_PL16     CRA     pre
lamich_PL4      CRA     pre
lamich_PL7      CRA     pre
lamich_PL9      CRA     pre
SM1083  FOR     post
SM1156  FOR     post
SM1204  FOR     post
SM1237  FOR     post
SM1270  FOR     post
SM1271  FOR     post
SM1272  FOR     post
SM1273  FOR     post
SRR2917289      FOR     pre
SRR2917290      FOR     pre
SRR2917291      FOR     pre
SRR2917292      FOR     pre
SRR2917293      FOR     pre
SRR2917294      FOR     pre
SRR2917295      FOR     pre
SRR2917296      FOR     pre
SRR2917297      FOR     pre
SRR2917298      FOR     pre
RHC097  PAR     post
RHC507  PAR     post
SM031   PAR     post
SM032   PAR     post
SM040   PAR     post
SM059   PAR     post
SM079   PAR     post
lamich_PARV1    PAR     pre
lamich_PARV2    PAR     pre
SRR2917329      PAR     pre
SRR2917330      PAR     pre
SRR2917331      PAR     pre
SRR2917332      PAR     pre
SRR2917333      PAR     pre
SRR2917334      PAR     pre
SRR2917335      PAR     pre
SRR2917336      PAR     pre
SRR2917337      PAR     pre
SRR2917338      PAR     pre





# plot it in R
cd  /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all

C <- as.matrix(read.table("pruned.cov")) # Reads estimated covariance matrix
# D <- as.matrix(read.table("output.selection")) # Reads PC based selection statistics
tab <- read.table("/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Treatment <- factor(labs$treatment)
labs$Sample <- factor(labs$sample)

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


# plot admixture
tbl=read.table("pruned.admix.3.Q")
pdf(file = "admix_pruned.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()




# plot PCA by each species
# CRA
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/cra

tail -n +2 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre_pops.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra.all.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post_pops.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra.all.pca.txt 

plink --bed /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra.all.pca.txt --out genolike_pruned_cra --make-bed 

pcangsd -p genolike_pruned_cra -o genolike_pruned_cra -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/cra

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.cra.txt
grep 'CRA' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.cra.txt

C <- as.matrix(read.table("genolike_pruned_cra.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.cra.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Treatment <- factor(labs$treatment)
labs$Sample <- factor(labs$sample)

pdf(file = "pca_pruned_cra_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_cra_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_cra_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("genolike_pruned_cra.admix.3.Q")
pdf(file = "admix_pruned_cra.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()


# FOR 

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/for

tail -n +2 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for.all.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post_pops.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for.all.pca.txt 

plink --bed /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for.all.pca.txt --out genolike_pruned_for --make-bed 

pcangsd -p genolike_pruned_for -o genolike_pruned_for -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/for

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.for.txt
grep 'FOR' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.for.txt

C <- as.matrix(read.table("genolike_pruned_for.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.for.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Treatment <- factor(labs$treatment)
labs$Sample <- factor(labs$sample)

pdf(file = "pca_pruned_for_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_for_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_for_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("genolike_pruned_for.admix.3.Q")
pdf(file = "admix_pruned_for.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()

# PAR

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/par

tail -n +2 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt | awk '{print $1,$1}' | awk '{print $1, $2, 0, 0, 0, -9}' > updated.fam

awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par.all.pca.txt 
awk 'BEGIN{OFS="\t"} {print $1, $1}'  /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post_pops.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par.all.pca.txt 

plink --bed /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bed --bim /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.bim --fam updated.fam --allow-extra-chr --snps-only 'just-acgt' --keep /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par.all.pca.txt --out genolike_pruned_par --make-bed 

pcangsd -p genolike_pruned_par -o genolike_pruned_par -t 12 -e 2 --selection --admix


# plot it in R
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/par

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.par.txt
grep 'PAR' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.txt >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.par.txt

C <- as.matrix(read.table("genolike_pruned_par.cov")) # Reads estimated covariance matrix
tab <- read.table("/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/sample_species_treatment.pca.par.txt", header = TRUE)
labs <- data.frame(tab)
# Plot PCA plot
e <- eigen(C)
PC1.PV.full = (e$values[1]/sum(e$values))*100
PC2.PV.full = (e$values[2]/sum(e$values))*100
PC1.PV = round(PC1.PV.full, digits = 2)
PC2.PV = round(PC2.PV.full, digits = 2)

labs$Species <- factor(labs$species)
labs$Treatment <- factor(labs$treatment)
labs$Sample <- factor(labs$sample)

pdf(file = "pca_pruned_par_species.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Species))
legend("topright", legend=levels(labs$Species), pch="o", col=1:nlevels(labs$Species))
dev.off()

pdf(file = "pca_pruned_par_treatment.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Treatment))
legend("topright", legend=levels(labs$Treatment), pch="o", col=1:nlevels(labs$Treatment))
dev.off()

pdf(file = "pca_pruned_par_sample.pdf", useDingbats=FALSE)
plot(e$vectors[,1:2], xlab=paste0("PC1 (Percent Variation =",PC1.PV,"%)"), ylab=paste0("PC2 (Percent Variation =",PC2.PV,"%)"), main="PCA", col=as.integer(labs$Sample))
legend("topright", legend=levels(labs$Sample), pch="o", col=1:nlevels(labs$Sample))
dev.off()

# plot admixture
tbl=read.table("genolike_pruned_par.admix.3.Q")
pdf(file = "admix_pruned_par.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 4)
barplot(t(as.matrix(tbl)), col=rainbow(3),
xlab="Individual #", ylab="Ancestry", border=NA)
dev.off()
