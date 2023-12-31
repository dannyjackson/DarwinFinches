#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide fst

args = commandArgs()

outDir <- args[6]
name <- args[7]

# for troubleshooting in interactive:
# outDir <- "/xdisk/mcnew/dannyjackson/finches/cra/fst"
# name <- "cra"

library(qqman)

fst<-read.table(paste0(outDir,"/",name,".chroms.windowed.weir.fst"), header=TRUE)

fstsubset<-fst[complete.cases(fst),]


xu <- mean(fstsubset$WEIGHTED_FST)
s <- sd(fstsubset$WEIGHTED_FST)
fstsubset$ZFST = (fstsubset$WEIGHTED_FST - xu)/s

SNP<-c(1: (nrow(fstsubset)))

mydf<-data.frame(SNP,fstsubset)


write.csv(mydf[ which(mydf$ZFST>='5'),], paste0(outDir,"/",name,".zfst_sig.csv"), row.names=FALSE)


pdf(file = paste0(outDir,"/",name,".chroms.fst.pdf"), width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()


pdf(file = paste0(outDir,"/",name,".chroms.zfst.pdf"), width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="ZFST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()