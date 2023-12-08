#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide dxy output from pixy

args = commandArgs()

outDir <- args[6]
name <- args[7]

# outDir <- "/xdisk/mcnew/dannyjackson/finches/cra/dxy"
# name <- "cra"

dxy.all <- read.table(paste0(outDir,"/cra.chroms.pixy_dxy.txt"),header=T)
dxy.subset<-dxy.all[complete.cases(dxy.all),]

SNP<-c(1: (nrow(dxy.subset)))

lower = min(dxy.subset$avg_dxy)
upper = max(dxy.subset$avg_dxy)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- dxy.subset$avg_dxy < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,dxy.subset)

pdf(file = paste0(outDir,"/",name,".dxy_hist.pdf"), width = 10, height = 5, useDingbats=FALSE)
hist(dxy.subset$avg_dxy,br=50)
dev.off()

pdf(file = paste0(outDir,"/",name,".dxy.pdf"), width = 20, height = 7, useDingbats=FALSE)

plot(avg_dxy ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], paste0(outDir,"/",name,".dxy_sig.csv"), row.names=FALSE)

