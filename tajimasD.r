#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide Tajima's D scan output from vcftools

args = commandArgs()

outDir <- args[6]
name <- args[7]



taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant


mydf<-data.frame(SNP,myBg,taj.subset)

sigdf <-  mydf[which(mydf$myBg),]

write.table(sigdf, file = "sigtd.tsv")

pdf(file = "pyrr_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "pyrr_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()