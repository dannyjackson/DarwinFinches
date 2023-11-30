#!/usr/bin/env Rscript
#Genome wide ABBA BABA test with block jackknife procedure

args = commandArgs()

outDir <- args[6]
name <- args[7]


R

pi.all <- read.table(paste0(outDir,"/",name,".chroms.windowed.pi"),header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = paste0(outDir,"/",name,".pre_pi_hist.pdf"), width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=50)
dev.off()

pdf(file = paste0(outDir,"/",name,".pre_pi.pdf"), width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

# this gives regions of significance 


write.csv(mydf[ which(mydf$myBg=='TRUE'),], paste0(outDir,"/",name,".pre_pi_sig.csv"), row.names=FALSE)

