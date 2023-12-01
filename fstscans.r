#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide fst

args = commandArgs()

outDir <- args[6]
name <- args[7]


pi.all <- read.table(paste0(outDir,"/fst/",name,".chroms.windowed.pi"),header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)




R
library(qqman)
fst<-read.table(paste0(outDir,"/fst/",name,"chroms.weir.formanhattan.fst"), header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = paste0(outDir,"/fst/",name,".chroms.fst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="SNP",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()
