#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide BayeScan

args = commandArgs()

outDir <- args[6]
name <- args[7]


library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR("/scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr.filtered.geno25.maf1.vcf")

pop_map <- read.table("/scratch/dnjacks4/cardinalis/to_parus/reference_lists/pyrr_bayescan_popfile.txt", header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn="pyrr.filtered.bsc")

mkdir filtered 

bayescan -n 5000 -burn 50000 -pr_odds 10 /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr/pyrr.filtered.bsc -od /scratch/dnjacks4/cardinalis/to_b10k/bayescan/pyrr/filtered/




