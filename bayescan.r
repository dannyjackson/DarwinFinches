#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide BayeScan

args = commandArgs()

outDir <- args[6]
name <- args[7]
pop1 <- args[8]
bayesmap <- args[9]

# for troubleshooting in interactive:
# outDir <- "/xdisk/mcnew/dannyjackson/finches/cra/bayescan"
# name <- "cra"

library(vcfR)
library(hierfstat)
library(adegenet)
library(ggplot2)
library(radiator)

vcf <- read.vcfR(paste0(pop1))

pop_map <- read.table(paste0(bayesmap), header=TRUE, stringsAsFactors = TRUE)

genind <- vcfR2genind(vcf)
genind@pop <- pop_map$STRATA
hierfstat <- genind2hierfstat(genind)


write.bayescan(hierfstat,fn=paste0(outDir,"/",name,".bsc"))