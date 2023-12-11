#!/usr/bin/env Rscript 
# Plotting and identifying significant regions of genome wide dxy output from pixy

args = commandArgs()

outDir <- args[6]
name <- args[7]

# outDir <- "/xdisk/mcnew/dannyjackson/finches/cra/pixy"
# name <- "cra"

## pi

pi.all <- read.table(paste0(outDir,"/pi/cra.chroms.pixy_pi.txt"),header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$avg_pi)
upper = max(pi.subset$avg_pi)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$avg_pi < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = paste0(outDir,"/pi/",name,".pi_hist.pdf"), width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$avg_pi,br=50)
dev.off()

pdf(file = paste0(outDir,"/pi/",name,".pi.pdf"), width = 20, height = 7, useDingbats=FALSE)

plot(avg_pi ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], paste0(outDir,"/pi/",name,".pi_sig.csv"), row.names=FALSE)


## fst

fst.all <- read.table(paste0(outDir,"/fst/cra.chroms.pixy_fst.txt"),header=T)
fst.subset<-fst.all[complete.cases(fst.all),]

SNP<-c(1: (nrow(fst.subset)))

lower = min(fst.subset$avg_wc_fst)
upper = max(fst.subset$avg_wc_fst)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- fst.subset$avg_wc_fst < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,fst.subset)

pdf(file = paste0(outDir,"/fst/",name,".fst_hist.pdf"), width = 10, height = 5, useDingbats=FALSE)
hist(fst.subset$avg_wc_fst,br=50)
dev.off()

pdf(file = paste0(outDir,"/fst/",name,".fst.pdf"), width = 20, height = 7, useDingbats=FALSE)

plot(avg_wc_fst ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$chromosome,
     las = 2)

dev.off()

# this gives regions of significance 

write.csv(mydf[ which(mydf$myBg=='TRUE'),], paste0(outDir,"/fst/",name,".fst_sig.csv"), row.names=FALSE)


## dxy 

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




