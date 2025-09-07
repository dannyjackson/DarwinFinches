samtools depth ${BAMDIR}/"$bird".realigned.bam -r NC_044601.1 >> ${OUTDIR}/datafiles/bamstats/"$bird"_depthstats.txt 

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/sexing_by_depth
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches/

BAMDIR=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment

while read -r bird; do
ls ${BAMDIR}/"$bird".realigned.bam >> ${OUTDIR}/referencelists/allsamplebams.txt 
done <  ${OUTDIR}/referencelists/allsamples.txt 

while read -r bird; do
echo $bird >> ${OUTDIR}/datafiles/bamstats/"$bird"_depthstats.txt 
samtools depth -r NC_044601.1 ${BAMDIR}/"$bird".realigned.bam >> "$bird"_depthstats.txt 
done <  ${OUTDIR}/referencelists/allsamples.txt 


while read -r bird; do
samtools idxstats ${BAMDIR}/"$bird".realigned.bam | grep 'NC' > "$bird".idxstats
done <  ${OUTDIR}/referencelists/allsamples.txt 


## source the R functions 
source("https://raw.githubusercontent.com/popgenDK/SATC/main/satcFunc.R")

## use the mclust library. Install if needed (install.package("mclust"))
library(mclust)

## read list of filenames
IDXFILE <- "idxindex.txt"
filenames <- scan(IDXFILE,what="sUp")

## read idx files
idx <- lapply(filenames,read.table,as.is=T)
names(idx) <- basename(filenames)

## Filter scafoolds (min 100kb ) and normalize using the M longest scaffold 
rFilt <- filterScaffold(dat=idx,minLength=1e5,M=5)

##plot normalized depth
pdf("plotDepth_allScaffs.pdf")
plotDepth(rFilt) 
dev.off()

## identify sex and sex scaffolds
sex <- sexDetermine(dat=rFilt, K=2, model="gaussian") 
x <- sexDetermine(dat=rFilt, K=2, model="gaussian") 


# Plot sex scaffolds' normalized depth per individual
pdf("plotSamples_sex.pdf")
plotSamples(sex)
dev.off()







# Plot clustering (sex inference)
pdf("plotGroup_sex.pdf")
plotGroup(sex)
dev.off()

pdf("plotScafs.pdf")
plotScafs(sex)
dev.off()

# Plot scaffold depths stratified by inferred sex
pdf("plotScafs_abnormal.pdf")
plotScafs(sex, abnormal = TRUE)
dev.off()


# Visualize Gaussian mixtures
pdf("plotUnc_gaussian.pdf")
plotUnc(sex)
dev.off()


## See the inferred status of each scaffold
 head(sex$SexScaffolds)
 
## See the inferred sex of each indiviual
cbind(names(idx),sex$sex)


[1,] "JP4481_all.idxstats"   "homomorphic"  
 [2,] "JP5410_all.idxstats"   "homomorphic"  
 [3,] "JP9655_all.idxstats"   "homomorphic"  
 [4,] "lamich_PARV1.idxstats" "homomorphic"  
 [6,] "lamich_PL15.idxstats"  "homomorphic"  
 [7,] "lamich_PL16.idxstats"  "homomorphic"  
 [8,] "lamich_PL4.idxstats"   "homomorphic"  
 [9,] "lamich_PL7.idxstats"   "homomorphic"  
[11,] "RHC097_all.idxstats"   "homomorphic"  
[12,] "RHC507_all.idxstats"   "homomorphic"  
[13,] "SM031_all.idxstats"    "homomorphic"  
[14,] "SM032_all.idxstats"    "homomorphic"  
[15,] "SM040_all.idxstats"    "homomorphic"  
[16,] "SM059_all.idxstats"    "homomorphic"  
[18,] "SM1067_all.idxstats"   "homomorphic"  
[19,] "SM1083.idxstats"       "homomorphic"  
[20,] "SM1156.idxstats"       "homomorphic"  
[21,] "SM1157_all.idxstats"   "homomorphic"  
[22,] "SM1200_all.idxstats"   "homomorphic"  
[23,] "SM1204_all.idxstats"   "homomorphic"  
[24,] "SM1231.idxstats"       "homomorphic"  
[25,] "SM1237.idxstats"       "homomorphic"  
[26,] "SM1240_all.idxstats"   "homomorphic"  
[27,] "SM1266_all.idxstats"   "homomorphic"  
[28,] "SM1270.idxstats"       "homomorphic"  
[29,] "SM1271.idxstats"       "homomorphic"  
[30,] "SM1272.idxstats"       "homomorphic"  
[31,] "SM1273.idxstats"       "homomorphic"  
[34,] "SRR2917291.idxstats"   "homomorphic"  
[37,] "SRR2917294.idxstats"   "homomorphic"  
[38,] "SRR2917295.idxstats"   "homomorphic"  
[40,] "SRR2917297.idxstats"   "homomorphic"  
[41,] "SRR2917298.idxstats"   "homomorphic"  
[42,] "SRR2917329.idxstats"   "homomorphic"  
[43,] "SRR2917330.idxstats"   "homomorphic"  
[44,] "SRR2917331.idxstats"   "homomorphic"  
[46,] "SRR2917333.idxstats"   "homomorphic"  
[47,] "SRR2917334.idxstats"   "homomorphic"  
[48,] "SRR2917335.idxstats"   "homomorphic"  
[50,] "SRR2917337.idxstats"   "homomorphic"  

 [5,] "lamich_PARV2.idxstats" "heteromorphic"
[10,] "lamich_PL9.idxstats"   "heteromorphic"
[17,] "SM079_all.idxstats"    "heteromorphic"
[32,] "SRR2917289.idxstats"   "heteromorphic"
[33,] "SRR2917290.idxstats"   "heteromorphic"
[35,] "SRR2917292.idxstats"   "heteromorphic"
[36,] "SRR2917293.idxstats"   "heteromorphic"
[39,] "SRR2917296.idxstats"   "heteromorphic"
[45,] "SRR2917332.idxstats"   "heteromorphic"
[49,] "SRR2917336.idxstats"   "heteromorphic"
[51,] "SRR2917338.idxstats"   "heteromorphic"

Males
SRR2917291
SRR2917294
SRR2917295
SRR2917297
SRR2917298
SRR2917329
SRR2917330
SRR2917331
SRR2917333
SRR2917334
SRR2917335
SRR2917337

Females
lamich_PARV2
lamich_PL9
SRR2917289
SRR2917290
SRR2917292
SRR2917293
SRR2917296
SRR2917332
SRR2917336
SRR2917338