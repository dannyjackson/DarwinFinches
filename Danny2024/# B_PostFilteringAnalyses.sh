# Post-Filtering analyses

/xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt

JP4481_all
JP5410_all
JP9655_all
SM1067_all
SM1157_all
SM1200_all
SM1231
SM1240_all
SM1266_all
lamich_PL15
lamich_PL16
lamich_PL4
lamich_PL7
lamich_PL9
SM1083
SM1156
SM1204_all
SM1237
SM1270
SM1271
SM1272
SM1273
SRR2917289
SRR2917290
SRR2917291
SRR2917292
SRR2917293
SRR2917294
SRR2917295
SRR2917296
SRR2917297
SRR2917298
RHC097_all
RHC507_all
SM031_all
SM032_all
SM040_all
SM059_all
SM079_all
lamich_PARV1
lamich_PARV2
SRR2917329
SRR2917330
SRR2917331
SRR2917332
SRR2917333
SRR2917334
SRR2917335
SRR2917336
SRR2917337
SRR2917338

/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_all.realigned.bam

sed -e 's/$/.realigned.bam/g' /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt | sed -e 's/^/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\//g' > /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt


#!/bin/bash

#SBATCH --job-name=angsdsaf
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/allsnps -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

sbatch angsd_saf.sh 
Submitted batch job 1962904


#!/bin/bash

#SBATCH --job-name=angsdsfs
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsfs.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/realSFS

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/allsnps.saf.idx -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -P 10 

sbatch realsfs.sh 
Submitted batch job 1963090

# 0.000000 0.000000 0.000000 2.098167 0.000000 45420.062922 899733.650117 403247.243583 397300.108820 344937.618533 289166.552720 295938.464004 231504.315574 211939.212108 207839.132847 181599.800006 169853.326312 163721.589632 150542.055493 142206.986045 135615.022760 126700.735216 118959.724721 118220.079971 118727.305085 103109.548810 93708.096243 104844.574559 227246.450338 230104.651333 135867.295839 96381.827280 88663.223918 90991.653759 86748.472458 77821.830350 73799.690522 71587.523039 69333.455410 64084.846319 62603.141033 61971.376724 60803.513211 56175.117846 55513.409503 55255.720983 56362.764873 56090.347777 43696.268285 36904.724551 40667.015591 136656.005797 13569.215890 24192.343766 46026.315194 56975.243720 48906.661691 44159.597999 45241.511999 46712.071655 48507.988607 50560.853237 47360.304327 48968.798590 49058.186792 48078.356951 42266.213452 41083.017891 37955.670343 35550.006393 34889.170062 33847.970957 32946.368834 30693.170356 30206.568736 29524.279589 28670.103056 28354.428025 28087.997856 27204.937591 26261.595754 26618.092407 25372.054628 25145.125818 24055.260330 23234.909041 24183.102904 23046.584115 22794.253175 22961.854370 21537.150068 22332.946396 21322.279568 21175.028686 21944.428536 19462.040127 30519.558343 66.729580 0.000000 0.000000 30.270714 154.659043 866.091847 




#!/bin/bash

#SBATCH --job-name=angsdsaf
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/postbams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/post/postsaf -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10


~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/prebams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/pre/presaf -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

sbatch angsd_saf_batch.sh 
Submitted batch job 1963119


# POST
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM059_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1083.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1156.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1231.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1237.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1270.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1271.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1272.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1273.realigned.bam


# PRE
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV2.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL15.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL16.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL4.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL7.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL9.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917289.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917290.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917291.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917292.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917293.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917294.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917295.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917296.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917297.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917298.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917329.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917330.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917331.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917332.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917333.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917334.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917335.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917336.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917337.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917338.realigned.bam

# The MAFs estimated in this step were used later to extract a list of private alleles in each batch of data (alleles with frequencies between 10% and 90% in one batch but smaller than 1% or larger than 99% in the other batch) (custom R script: https://github.com/therkildsen-lab/batch-effect/blob/main/markdown/degradation.md#extract-private-alleles-and-examine-proportion-of-base-substitutions).


library(tidyverse)
library(ggplot)

maf_pre <- read_tsv("/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/pre/presaf.mafs.gz") %>%
  transmute(lg = chromo, position = position, major=major, minor = minor, pre_maf = knownEM, pre_nind=nInd)
maf_post <- read_tsv("/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/post/postsaf.mafs.gz")%>%
  transmute(lg = chromo, position = position, major=major, minor = minor, post_maf = knownEM, post_nind=nInd)

maf_joined <- inner_join(maf_pre, maf_post) %>%
  mutate(delta = abs(pre_maf- post_maf))



private_alleles <- bind_rows((filter(maf_joined, pre_maf<0.01 | pre_maf>0.99) %>% filter(post_maf>0.1 & post_maf<0.9) %>% transmute(major=major, minor=minor, batch = "post samples")),
          (filter(maf_joined, post_maf<0.01 | post_maf>0.99) %>% filter(pre_maf>0.1 & pre_maf<0.9) %>% transmute(major=major, minor=minor, batch = "pre samples"))) %>%
  mutate(base_substitution = str_c(major, "-to-", minor)) %>%
  mutate(base_substitution = case_when(
    base_substitution %in% c("A-to-C", "T-to-G") ~ "A-to-C\nT-to-G", 
    base_substitution %in% c("A-to-G", "T-to-C") ~ "A-to-G\nT-to-C", 
    base_substitution %in% c("A-to-T", "T-to-A") ~ "A-to-T\nT-to-A", 
    base_substitution %in% c("C-to-A", "G-to-T") ~ "C-to-A\nG-to-T", 
    base_substitution %in% c("C-to-G", "G-to-C") ~ "C-to-G\nG-to-C", 
    base_substitution %in% c("C-to-T", "G-to-A") ~ "C-to-T\nG-to-A"
  )) 

p_b <- private_alleles %>%
  group_by(base_substitution, batch) %>% 
  count() %>% 
  ungroup() %>% 
  group_by(batch) %>%
  mutate(frequency = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=base_substitution, y=frequency, fill=batch, group=batch)) +
  #geom_line() + 
  #geom_point() +
  geom_col(position = "dodge", color="black")+
  scale_fill_viridis_d(begin = 0.25, end=0.75) +
  ylim(c(NA, 0.4)) +
  labs(x="type of base substitution", y="frequency") +
  theme(legend.position = "top",
        panel.grid =element_blank(),
        axis.line = element_line())


pdf(file = "privatealleles_1.pdf", width = 15, height =10, useDingbats=FALSE)
print(p_b)
dev.off()

# estimate individual heterozygosity with and without transitions
/xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt


#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j

while read -r finch;
do
echo "$finch" >> /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/progress.txt

~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/"$finch"

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt 

sbatch angsdsaf_indv.sh 
Submitted batch job 9910629

#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=3:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j

~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917337.realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/SRR2917337

~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917338.realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/SRR2917338

sbatch angsdsaf_indv2.sh 
Submitted batch job 9919022

#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j

while read -r finch;
do
echo "$finch" >> /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/progress.txt

~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/"$finch" -noTrans 1

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt 

sbatch angsdsaf_indv_notrans.sh 
Submitted batch job 9919414

#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=3:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j



~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917337.realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/SRR2917337 -noTrans 1

~/programs/angsd/angsd -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917338.realigned.bam -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -dosaf 1 -gl 1 -setMinDepthInd 4 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -nThreads 10 -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/SRR2917338 -noTrans 1

sbatch angsdsaf_indv_notrans2.sh 
Submitted batch job 9926857




#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j


while read -r finch;
do
echo "$finch" >> /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/progress.txt

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/"$finch".saf.idx > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/allsites/"$finch".est.ml

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt 


sbatch angsdsfs_indv.sh 
Submitted batch job 1973271

#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j


while read -r finch;
do

~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/"$finch".saf.idx > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/heterozygosity/no_transitions/"$finch".est.ml

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt 

sbatch angsdsfs_indv2.sh 
Submitted batch job 1973272

#!/usr/bin/env Rscript 
# Calculating heterozygosity from realSFS output

args = commandArgs()

finch <- args[6]

a<-scan(paste0(finch,".est.ml"), quiet = TRUE)

het = a[2]/sum(a)

cat(finch, het, "\n")



while read -r finch;
do

Rscript ~/programs/Rscripts/het.r $finch

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamples.txt 

# FST calculation -- not using LD filter
# Calculate saf files for each population
# CRA all
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1231.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL15.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL16.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL4.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL7.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL9.realigned.bam

#!/bin/bash

#SBATCH --job-name=angsdsaf_indv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_indv.%j

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_all_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/all/all -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

sbatch cra_all.sh 
Submitted batch job 1973604

# CRA pre
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL15.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL16.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL4.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL7.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PL9.realigned.bam


# CRA post
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1231.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266.realigned.bam



# PAR all
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM059_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV2.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917329.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917330.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917331.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917332.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917333.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917334.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917335.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917336.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917337.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917338.realigned.bam

# PAR pre
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/lamich_PARV2.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917329.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917330.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917331.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917332.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917333.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917334.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917335.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917336.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917337.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917338.realigned.bam

# PAR post
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM059_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_all.realigned.bam

# FOR all
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1083.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1156.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1237.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1270.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1271.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1272.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1273.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917289.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917290.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917291.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917292.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917293.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917294.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917295.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917296.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917297.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917298.realigned.bam

# FOR pre
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917289.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917290.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917291.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917292.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917293.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917294.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917295.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917296.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917297.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SRR2917298.realigned.bam

# FOR post
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1083.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1156.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1237.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1270.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1271.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1272.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1273.realigned.bam

#!/bin/bash

#SBATCH --job-name=angsdsaf_pops_forpost
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=25:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdsaf_pops_forpost.%j

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_post -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/par_all_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/all/par_all -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/par_pre_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/par_post_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_all_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/all/for_all -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pre_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 10

sbatch angsdsaf_pops.sh 
Submitted batch job 9931679



# FOR 

#this is with 2pops
#first calculate per pop saf for each population (done above)
/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx
#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx > pre.post.ml
#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx -sfs pre.post.ml -fstout here
#get the global estimate
~/programs/angsd/misc/realSFS fst stats here.fst.idx 



	-> FST.Unweight[nObs:8400488]:0.030847 Fst.Weight:0.031077

#below is not tested that much, but seems to work


~/programs/angsd/misc/realSFS fst stats2 here.fst.idx  -win 50000 -step 50000 >slidingwindow

~/programs/angsd/misc/realSFS fst stats2 here.fst.idx -win 1 -step 1 >slidingwindow_singlesnps

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps >> slidingwindow_singlesnps_fst.txt

sed -i 's/NC_//g' slidingwindow_singlesnps_fst.txt 


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow >> slidingwindow_fst.txt

sed -i 's/NC_//g' slidingwindow_fst.txt 
# sed -i 's/NW_//g' slidingwindow_fst.txt 
# sed -i 's/.1//g' slidingwindow_fst.txt 

# just some quick looks at the output:
fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
max(fst$fst)
subset(fst, fst > 0.7)
                                                       region         chr
3936651 (187634,187634)(25267982,25267982)(25267982,25267983) NC_044575.1
5445766 (134702,134702)(16843164,16843164)(16843164,16843165) NC_044580.1
5445767 (134703,134703)(16843166,16843166)(16843166,16843167) NC_044580.1
8139248 (279030,279030)(48092710,48092710)(48092710,48092711) NC_044601.1
          midPos Nsites      fst
3936651 25267982      2 0.824873
5445766 16843164      2 0.716541
5445767 16843166      2 0.744155
8139248 48092710      2 0.701127

# Only one region showed up twice, NC_044580.1 positions 16843164 and 16843166. Both of these are within the gene: synemin SYNM (NC_044580.1:16846280-16866830)
This is an intermediate filament protein originally found in avian smooth muscles. Gene Cards says it is a "cytoskeletal protein that confers resistance to mechanical stress"

library(qqman)
fst <- read.csv('slidingwindow_fst.txt', sep ='\t')

pdf(file = "manhattan_windowed_for.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()

library(qqman)

fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')

png(file = "manhattan_snps_par.png", width = 1425, height =975)
manhattan(fst, chr="chr", bp="midPos", snp="Nsites", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()

# CRA
#this is with 2pops
#first calculate per pop saf for each population (done above)
/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx
#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx > pre.post.ml
#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx -sfs pre.post.ml -fstout here

#get the global estimate
~/programs/angsd/misc/realSFS fst stats here.fst.idx 

# cra:
-> FST.Unweight[nObs:8281787]:0.024960 Fst.Weight:0.036349

# for:
	-> FST.Unweight[nObs:8400488]:0.030847 Fst.Weight:0.031077

#below is not tested that much, but seems to work

~/programs/angsd/misc/realSFS fst stats2 here.fst.idx  -win 50000 -step 10000 >slidingwindow
~/programs/angsd/misc/realSFS fst stats2 here.fst.idx -win 1 -step 1 >slidingwindow_singlesnps

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps >> slidingwindow_singlesnps_fst.txt

sed -i 's/NC_//g' slidingwindow_singlesnps_fst.txt 


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow >> slidingwindow_fst.txt

sed -i 's/NC_//g' slidingwindow_fst.txt 
# sed -i 's/NW_//g' slidingwindow_fst.txt 
# sed -i 's/.1//g' slidingwindow_fst.txt 

# just some quick looks at the output:
fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
max(fst$fst)
subset(fst, fst > 0.7)

                                                           region     chr
335661      (335661,335661)(37907290,37907290)(37907290,37907291) 44571.1
335663      (335663,335663)(37907336,37907336)(37907336,37907337) 44571.1
3048343 (820574,820574)(104906414,104906414)(104906414,104906415) 44573.1
3048345 (820576,820576)(104906573,104906573)(104906573,104906574) 44573.1
3048362 (820593,820593)(104908569,104908569)(104908569,104908570) 44573.1
3048369 (820600,820600)(104909982,104909982)(104909982,104909983) 44573.1
3048373 (820604,820604)(104910683,104910683)(104910683,104910684) 44573.1
4083631     (356727,356727)(50553761,50553761)(50553761,50553762) 44575.1
4416458     (225338,225338)(27896405,27896405)(27896405,27896406) 44576.1
6735999     (548241,548241)(65558014,65558014)(65558014,65558015) 44586.1
6842675           (58407,58407)(6978972,6978972)(6978972,6978973) 44587.1
7861887     (109171,109171)(13085797,13085797)(13085797,13085798) 44601.1
           midPos Nsites      fst
335661   37907290      2 0.724117
335663   37907336      2 0.715765
3048343 104906414      2 0.769188
3048345 104906573      2 0.738477
3048362 104908569      2 0.729066
3048369 104909982      2 0.865019
3048373 104910683      2 0.728000
4083631  50553761      2 0.770351
4416458  27896405      2 0.720562
6735999  65558014      2 0.700035
6842675   6978972      2 0.716074
7861887  13085797      2 0.789113

# NC_044571.1   pos: 37907290 and 37907336
# this is not within a gene, but is just after:
NC_044571.1:37901151-37906295 testis-expressed protein 30 isoform X1	TEX30
# NCBI says: "Predicted to enable hydrolase activity."
# Gene cards says: "Diseases associated with TEX30 include Spermatogenic Failure 25 and Hypotonia-Cystinuria Syndrome"

# and before: 
NC_044571.1:37937016-37945300 LOW QUALITY PROTEIN: protein-lysine methyltransferase METTL21C	METTL21C
# NCBI says: "Enables heat shock protein binding activity and protein-lysine N-methyltransferase activity. Involved in protein methylation. Located in nucleus. Part of protein-containing complex."
# Gene cards says: "Diseases associated with METTL21C include Myostatin-Related Muscle Hypertrophy and Severe Congenital Neutropenia 6."

# 44573.1       pos: 104906414, 104906573, 104908569, 104909982, 104910683
Between these two:
NC_044573.1:104387453-104741091 glutamate receptor ionotropic, kainate 2	GRIK2
# NCBI says: "Glutamate receptors are the predominant excitatory neurotransmitter receptors in the mammalian brain and are activated in a variety of normal neurophysiologic processes. "
# gene cards says: "Diseases associated with GRIK2 include Neurodevelopmental Disorder With Impaired Language And Ataxia And With Or Without Seizures and Intellectual Developmental Disorder, Autosomal Recessive 6. Among its related pathways are Presynaptic function of Kainate receptors and Transmission across Chemical Synapses."
NC_044573.1:105878099-105927569 E3 ubiquitin-protein ligase HACE1 isoform X1	HACE1
# NCBI says: "The encoded protein is involved in specific tagging of target proteins, leading to their subcellular localization or proteasomal degradation. The protein is a potential tumor suppressor and is involved in the pathophysiology of several tumors, including Wilm's tumor."
# Gene cards says: "Diseases associated with HACE1 include Spastic Paraplegia And Psychomotor Retardation With Or Without Seizures and Charcot-Marie-Tooth Disease, Axonal, Type 2W. Among its related pathways are Class I MHC mediated antigen processing and presentation and Innate Immune System."

module load R


library(qqman)
fst <- read.csv('slidingwindow_fst.txt', sep ='\t')

pdf(file = "manhattan.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()


library(qqman)
fst <- read.csv('slidingwindow_fst.txt', sep ='\t')

pdf(file = "manhattan_windowed_cra.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()











# PAR

#this is with 2pops
#first calculate per pop saf for each population (done above)
/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx
#calculate the 2dsfs prior
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx > pre.post.ml
#prepare the fst for easy window analysis etc
~/programs/angsd/misc/realSFS fst index /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx -sfs pre.post.ml -fstout here
#get the global estimate
~/programs/angsd/misc/realSFS fst stats here.fst.idx 

# par: 
-> FST.Unweight[nObs:8527720]:0.028750 Fst.Weight:0.028649

# cra:
-> FST.Unweight[nObs:8281787]:0.024960 Fst.Weight:0.036349

# for:
	-> FST.Unweight[nObs:8400488]:0.030847 Fst.Weight:0.031077


#below is not tested that much, but seems to work

~/programs/angsd/misc/realSFS fst stats2 here.fst.idx  -win 50000 -step 10000 >slidingwindow
~/programs/angsd/misc/realSFS fst stats2 here.fst.idx -win 1 -step 1 >slidingwindow_singlesnps

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps >> slidingwindow_singlesnps_fst.txt

sed -i 's/NC_//g' slidingwindow_singlesnps_fst.txt 


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow >> slidingwindow_fst.txt

sed -i 's/NC_//g' slidingwindow_fst.txt 
# sed -i 's/NW_//g' slidingwindow_fst.txt 
# sed -i 's/.1//g' slidingwindow_fst.txt 

# just some quick looks at the output:
fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
max(fst$fst)
subset(fst, fst > 0.7)

                                                       region     chr   midPos
5125132 (229970,229970)(27663929,27663929)(27663929,27663930) 44578.1 27663929
5125133 (229971,229971)(27663934,27663934)(27663934,27663935) 44578.1 27663934
        Nsites      fst
5125132      2 0.710722
5125133      2 0.701462

subset(fst, fst > 0.7)

# NC_044578.1:27663929-27663934
Between the following two genes:
NC_044578.1:27643755-27644801 LOW QUALITY PROTEIN: immediate early response gene 5 protein	IER5
# NCBI says: "This gene encodes a protein that is similar to other immediate early response proteins. In the mouse, a similar gene may play an important role in mediating the cellular response to mitogenic signals. Studies in rats found the expression of a similar gene to be increased after waking and sleep deprivation."
# Gene cards doesn't say anything useful

NC_044578.1:27705862-27829414 voltage-dependent R-type calcium channel subunit alpha-1E isoform X1	CACNA1E
# NCBI says: "Voltage-dependent calcium channels are multisubunit complexes consisting of alpha-1, alpha-2, beta, and delta subunits in a 1:1:1:1 ratio. These channels mediate the entry of calcium ions into excitable cells, and are also involved in a variety of calcium-dependent processes, including muscle contraction, hormone or neurotransmitter release, gene expression, cell motility, cell division and cell death. This gene encodes the alpha-1E subunit of the R-type calcium channels, which belong to the 'high-voltage activated' group that maybe involved in the modulation of firing patterns of neurons important for information processing."
# Gene cards says "Diseases associated with CACNA1E include Developmental And Epileptic Encephalopathy 69 and Van Der Woude Syndrome 1. Among its related pathways are DREAM Repression and Dynorphin Expression and TCR Signaling (Qiagen)."



> subset(fst, fst > 0.65)

                                                           region     chr
757361      (757361,757361)(87171218,87171218)(87171218,87171219) 44571.1
757362      (757362,757362)(87171219,87171219)(87171219,87171220) 44571.1
2568818     (304218,304218)(38140696,38140696)(38140696,38140697) 44573.1
3046746     (782146,782146)(98996648,98996648)(98996648,98996649) 44573.1
3128069 (863469,863469)(108909526,108909526)(108909526,108909527) 44573.1
3979790     (189914,189914)(25219440,25219440)(25219440,25219441) 44575.1
4494346     (229325,229325)(27901632,27901632)(27901632,27901633) 44576.1
5125132     (229970,229970)(27663929,27663929)(27663929,27663930) 44578.1
5125133     (229971,229971)(27663934,27663934)(27663934,27663935) 44578.1
5349271     (196554,196554)(23290519,23290519)(23290519,23290520) 44579.1
           midPos Nsites      fst
757361   87171218      2 0.667568
757362   87171219      2 0.668538
2568818  38140696      2 0.672893
3046746  98996648      2 0.692131
3128069 108909526      2 0.662883
3979790  25219440      2 0.679678
4494346  27901632      2 0.675304
5125132  27663929      2 0.710722
5125133  27663934      2 0.701462
5349271  23290519      2 0.678239
# additional region of:
# NC_044571.1:87171218-87171219
# between these two genes:
NC_044571.1:87095605-87147397 transcriptional enhancer factor TEF-3	TEAD4
# NCBI says: "This gene product is a member of the transcriptional enhancer factor (TEF) family of transcription factors, which contain the TEA/ATTS DNA-binding domain. It is preferentially expressed in the skeletal muscle, and binds to the M-CAT regulatory element found in promoters of muscle-specific genes to direct their gene expression."
# Gene cards says: "Diseases associated with TEAD4 include Sveinsson Chorioretinal Atrophy. Among its related pathways are Gene expression (Transcription) and ERK Signaling."

NC_044571.1:87153063-87322660 tetraspanin-9 isoform X1	TSPAN9
# NCBI says: "The protein encoded by this gene is a member of the transmembrane 4 superfamily, also known as the tetraspanin family. Most of these members are cell-surface proteins that are characterized by the presence of four hydrophobic domains. The proteins mediate signal transduction events that play a role in the regulation of cell development, activation, growth and motility. "
# Gene cards didn't say anything useful

library(qqman)
fst <- read.csv('slidingwindow_fst.txt', sep ='\t')

pdf(file = "manhattan_windowed_par.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()




# make manhattan plot of relevant chroms:

echo -e 'region\tchr\tmidPos\tNsites\tfst' > relevantchroms_singlesnps.txt
grep 'NC_044578.1' slidingwindow_singlesnps >> relevantchroms_singlesnps.txt
grep 'NC_044571.1' slidingwindow_singlesnps >> relevantchroms_singlesnps.txt
sed -i 's/NC_//g' relevantchroms_singlesnps.txt 


library(qqman)
fst <- read.csv('relevantchroms_singlesnps.txt', sep ='\t')

png(file = "manhattan_relevantchroms_snps_par.png", width = 1425, height =975)
manhattan(fst, chr="chr", bp="midPos", snp="Nsites", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()



# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('relevantchroms_singlesnps.txt', sep ='\t')

df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.relevant.png", width=1425, height=975)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=2) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  # scale_y_continuous(expand = c(0, 0), limits = ylims) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle("title") +
  labs(x = "Chromosome") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.7)
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()

#!/bin/bash

#SBATCH --job-name=for.fstplot
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=25:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.for.fstplot.%j

#!/usr/bin/env Rscript
#Fst.R
module load R

Rscript /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/fst.rscript.r

sbatch fstplot.sh 
Submitted batch job 1978118

library(qqman)
fst <- read.csv('slidingwindow_fst.txt', sep ='\t')

pdf(file = "manhattan_allsnps.pdf", width = 20, height =7, useDingbats=FALSE)
manhattan(fst, chr="chr", bp="Nsites", snp="midPos", p="fst", logp=FALSE, ylab = "FST", cex = 0.5)
dev.off()








# FST PAR


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps | grep -v 'NC_044601.1' >> slidingwindow_singlesnps_fst.txt

echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow | grep -v 'NC_044601.1' >> slidingwindow_fst.txt

sed -i 's/NC_044571.1/1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044601.1/Z/g' slidingwindow_singlesnps_fst.txt


sed -i 's/NC_044571.1/1/g' slidingwindow_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_fst.txt
sed -i 's/NC_044601.1/Z/g' slidingwindow_fst.txt


interactive -a mcnew -m 50

#!/usr/bin/env Rscript 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.fst.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


#!/bin/bash

#SBATCH --job-name=par.fstplot
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=1:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.par.fstplot.%j

module load R

Rscript /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/fstplot.R

sbatch par.fstplot.sh 
Submitted batch job 1980722


# sliding window plot

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_fst.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.fst.windowed.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  # geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()



# FST FOR


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps | grep -v 'NC_044601.1' >> slidingwindow_singlesnps_fst.txt


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow | grep -v 'NC_044601.1' >> slidingwindow_fst.txt


sed -i 's/NC_044571.1/1/g' slidingwindow_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_fst.txt
sed -i 's/NC_044601.1/Z/g' slidingwindow_fst.txt




#!/usr/bin/env Rscript 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.fst.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


#!/bin/bash

#SBATCH --job-name=for.fstplot
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=1:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.for.fstplot.%j

module load R

Rscript /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/fstplot.R

sbatch for.fstplot.sh 
Submitted batch job 1980724






# windowed



#!/usr/bin/env Rscript 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_fst.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.fst.windowed.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  # geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# FST CRA


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_singlesnps_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow_singlesnps | grep -v 'NC_044601.1' >> slidingwindow_singlesnps_fst.txt


echo -e 'region\tchr\tmidPos\tNsites\tfst' > slidingwindow_fst.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' slidingwindow | grep -v 'NC_044601.1' >> slidingwindow_fst.txt

sed -i 's/NC_044571.1/1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_singlesnps_fst.txt
sed -i 's/NC_044601.1/Z/g' slidingwindow_singlesnps_fst.txt


sed -i 's/NC_044571.1/1/g' slidingwindow_fst.txt
sed -i 's/NC_044572.1/2/g' slidingwindow_fst.txt
sed -i 's/NC_044573.1/3/g' slidingwindow_fst.txt
sed -i 's/NC_044574.1/4/g' slidingwindow_fst.txt
sed -i 's/NC_044575.1/5/g' slidingwindow_fst.txt
sed -i 's/NC_044576.1/6/g' slidingwindow_fst.txt
sed -i 's/NC_044577.1/7/g' slidingwindow_fst.txt
sed -i 's/NC_044578.1/8/g' slidingwindow_fst.txt
sed -i 's/NC_044579.1/9/g' slidingwindow_fst.txt
sed -i 's/NC_044580.1/10/g' slidingwindow_fst.txt
sed -i 's/NC_044581.1/11/g' slidingwindow_fst.txt
sed -i 's/NC_044582.1/12/g' slidingwindow_fst.txt
sed -i 's/NC_044583.1/13/g' slidingwindow_fst.txt
sed -i 's/NC_044584.1/14/g' slidingwindow_fst.txt
sed -i 's/NC_044585.1/15/g' slidingwindow_fst.txt
sed -i 's/NC_044586.1/1.1/g' slidingwindow_fst.txt
sed -i 's/NC_044587.1/17/g' slidingwindow_fst.txt
sed -i 's/NC_044588.1/18/g' slidingwindow_fst.txt
sed -i 's/NC_044589.1/19/g' slidingwindow_fst.txt
sed -i 's/NC_044590.1/20/g' slidingwindow_fst.txt
sed -i 's/NC_044591.1/21/g' slidingwindow_fst.txt
sed -i 's/NC_044592.1/22/g' slidingwindow_fst.txt
sed -i 's/NC_044593.1/23/g' slidingwindow_fst.txt
sed -i 's/NC_044594.1/24/g' slidingwindow_fst.txt
sed -i 's/NC_044595.1/25/g' slidingwindow_fst.txt
sed -i 's/NC_044596.1/26/g' slidingwindow_fst.txt
sed -i 's/NC_044597.1/27/g' slidingwindow_fst.txt
sed -i 's/NC_044598.1/28/g' slidingwindow_fst.txt
sed -i 's/NC_044599.1/29/g' slidingwindow_fst.txt
sed -i 's/NC_044600.1/4.1/g' slidingwindow_fst.txt
sed -i 's/NC_044601.1/Z/g' slidingwindow_fst.txt



#!/usr/bin/env Rscript 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.fst.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  # geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


#!/bin/bash

#SBATCH --job-name=cra.fstplot
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=1:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.for.fstplot.%j

module load R

Rscript /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/fstplot.R

sbatch cra.fstplot.sh 
Submitted batch job 2043854

# windowed



#!/usr/bin/env Rscript 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.fst.windowed.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  # geom_hline(yintercept = 0.7) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()




# July 10th
# Trying to remember what the process is for determining a cutoff with these various populations statistics. I believe it has to do with SE but I need to find the paper that I was going to base my analyses off.

# From this paper:
# Novel signals of adaptive genetic variation in northwestern
# Clucas et al 2019
# We identified the genomic windows that were in the upper 99.9th percentile of the windowed FST distribution for each comparison, after excluding the inversions, and extracted all gene annotations from the gadMor2 reference genome (using the filtered gene set, which includes only putatively reliable annotations (Tørresen et al., 2017)) that were within 15 kb of the center of each window, thus in- vestigating a 30 kb window in total.
# I think I could use the 99th percentile logic for all analyses here too...


# subsetting top 1% of snps
subset(fst, fst > 0.7)

interactive -a mcnew -m 50

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
fst <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
max(fst$fst)
min(fst$fst)

ordered_fst <- fst %>% 
 # desc orders from largest to smallest
 arrange(desc(fst)) 

# 7752716 snps, so top 0.1% would be 7753

outlier_fst_disorder <- ordered_fst[1:7753,]

outlier_fst <- outlier_fst_disorder %>% arrange(chr, midPos)

min(outlier_fst_disorder$fst)
# 0.362058
max(outlier_fst_disorder$fst)
# 0.865019
write.csv(outlier_fst, "outlierfst.csv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_singlesnps_fst.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(midPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, midPos) %>%
  mutate( BPcum=midPos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.fst.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(fst))) +
  # Show all points
  geom_point(aes(color=as.factor(chr)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "FST") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.362058)
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midPos), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()






# moved tajima's d and fst codes to separate documents to keep better organized
# nucleotide diversity is next!


~/programs/angsd/misc/realSFS print pre/for_pre.saf.idx post/for_post.saf.idx | cut -f 1-2 > intersect.prepost.txt

~/programs/angsd/angsd sites index intersect.prepost.txt

NSITES=`wc -l intersect.prepost.txt | cut -f 1 -d " "`
echo $NSITES
# 8400488
zcat pre/for_pre.saf.gz > pre/for_pre.saf
zcat post/for_post.saf.gz > post/for_post.saf

zcat Results/TSI.saf.gz > Results/TSI.saf
zcat Results/PEL.saf.gz > Results/PEL.saf
NSITES=`wc -l intersect.prepost.txt | cut -f 1 -d " "` # if not already done
$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles Results/TSI.saf Results/PEL.saf -nsites $NSITES -nind 10 10 -outfile Results/TSI.PEL.stats.txt


#!/bin/bash

#SBATCH --job-name=dxy
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dxy.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/


~/programs/angsd/misc/realSFS print /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.idx | cut -f 1-2 > intersect.prepost.txt

~/programs/angsd/angsd sites index intersect.prepost.txt
NSITES=`wc -l intersect.prepost.txt | cut -f 1 -d " "`
echo $NSITES
# 8400488



# zcat Results/TSI.saf.gz > Results/TSI.saf
# zcat Results/PEL.saf.gz > Results/PEL.saf

$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles Results/TSI.saf Results/PEL.saf -nsites $NSITES -nind 10 10 -outfile Results/TSI.PEL.stats.txt

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.dxy.saf

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.dxy.saf

~/programs/ngsTools/ngsPopGen/ngsStat -npop 2 -postfiles /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.saf -nsites $NSITES -nind 8 10 -outfile /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre_post_dxy.txt

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.mafs

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.mafs


#!/bin/bash

#SBATCH --job-name=dxy.fortis
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=2:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dxy.fortis.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/dxy

module load R

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.mafs 

sbatch for.dxy.sh 
Submitted batch job 2045468

-r NC_044571.1


-t 8400488

# July 15th

zcat pre/presaf.saf.gz > pre/presaf.saf
zcat post/postsaf.saf.gz > post/postsaf.saf
NSITES=`wc -l Data/intersect.txt | cut -f 1 -d " "` # if not already done
$NGSTOOLS/ngsPopGen/ngsStat -npop 2 -postfiles Results/TSI.saf Results/PEL.saf -nsites $NSITES -nind 10 10 -outfile Results/TSI.PEL.stats.txt


sbatch forpost_angsd.sh 
Submitted batch job 2045372



# July 16th

/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/dxy
-r NC_044571.1

#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy


~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs 
# -r NC_044571.1

mv angsdput.mafs.gz for_post.mafs.gz

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pre_bams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs 
# -r NC_044571.1

mv angsdput.mafs.gz for_pre.mafs.gz


/xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/for_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/for_post.mafs


Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/for_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/for_post.mafs -t 604661

[1] "Global dxy is: 232582.909011818"
[1] "Global per site Dxy is: 0.38465009155844"


# try with previously made  mafs files
#!/bin/bash

#SBATCH --job-name=dxy_for
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dxy_for.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/for

module load R

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/pre/for_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.mafs -t 8400488

sbatch dxy_for.sh 
Submitted batch job 10116758

#!/bin/bash

#SBATCH --job-name=dxy_par
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dxy_par.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/par

module load R

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.mafs -t 8527720

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/pre/par_pre.mafs

zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.mafs

sbatch dxy_par.sh 
Submitted batch job 10116756

#!/bin/bash

#SBATCH --job-name=dxy_cra
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dxy_cra.%j

cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/dxy/cra

module load R

Rscript ~/programs/ngsTools/ngsPopGen/scripts/calcDxy.R -p /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.mafs -q /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.mafs -t 8281787


sbatch dxy_cra.sh 
Submitted batch job 10116664



zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/pre/cra_pre.mafs
zcat /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.mafs.gz > /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.mafs





# Plot dxy 

# cra
echo -e 'chromo\tposition\tdxy' > slidingwindow_dxy.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' Dxy_persite.txt >> slidingwindow_dxy.txt


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('slidingwindow_dxy.txt', sep ='\t')
max(dxy$dxy)
min(dxy$dxy)

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(dxy)) 

nrow(dxy)
# 8185302 snps, so top 0.1% would be 8185

outlier_dxy_disorder <- ordered_dxy[1:8185,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromo, position)

min(outlier_dxy_disorder$dxy)
# 0.6599634
max(outlier_dxy_disorder$dxy)
# 1
write.csv(outlier_dxy, "outlierdxy.csv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.dxy.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.6599634) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# for 

echo -e 'chromo\tposition\tdxy' > slidingwindow_dxy.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' Dxy_persite.txt >> slidingwindow_dxy.txt


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('slidingwindow_dxy.txt', sep ='\t')
max(dxy$dxy)
min(dxy$dxy)

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(dxy)) 

nrow(dxy)
# 8291664 snps, so top 0.1% would be 8291

outlier_dxy_disorder <- ordered_dxy[1:8291,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromo, position)

min(outlier_dxy_disorder$dxy)
# 0.6633972
max(outlier_dxy_disorder$dxy)
# 1
write.csv(outlier_dxy, "outlierdxy.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.dxy.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.6633972) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# par 

echo -e 'chromo\tposition\tdxy' > slidingwindow_dxy.txt
#tail -n+2 slidingwindow >> slidingwindow_fst.txt 
grep 'NC_' Dxy_persite.txt >> slidingwindow_dxy.txt


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
dxy <- read.csv('slidingwindow_dxy.txt', sep ='\t')
max(dxy$dxy)
min(dxy$dxy)

ordered_dxy <- dxy %>% 
 # desc orders from largest to smallest
 arrange(desc(dxy)) 

nrow(dxy)
# 8291664 snps, so top 0.1% would be 8291

outlier_dxy_disorder <- ordered_dxy[1:8291,]

outlier_dxy <- outlier_dxy_disorder %>% arrange(chromo, position)

min(outlier_dxy_disorder$dxy)
# 0.6111659
max(outlier_dxy_disorder$dxy)
# 1
write.csv(outlier_dxy, "outlierdxy.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chromo) %>% 
  summarise(chr_len=max(position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chromo"="chromo")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chromo, position) %>%
  mutate( BPcum=position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chromo) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.dxy.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(dxy))) +
  # Show all points
  geom_point(aes(color=as.factor(chromo)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromo, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "dxy") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 0.6633972) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


# pixy autosomes finished, check in on those


# SweeD
I need to make one big file with col1 as chromosome

# zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_urban.txt


grep 'Position' SweeD_Report.cra_pre_NC_044571.1 | head -n 1 > cra_pre.SweeD_Report.header


grep 'Position' SweeD_Report.cra_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.cra_pre

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.cra_pre
    fi
done < SweeD_Report.cra_pre_"$chrom" 
done < chroms.txt

grep 'Position' SweeD_Report.urban_noca | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_noca_urban.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.noca_urban_manhattan_scaffnames
  fi
done < SweeD_Report.noca_urban_manhattan

sed -i 's/VYXE//g' SweeD_Report.noca_urban_manhattan_scaffnames





#!/bin/bash

#SBATCH --job-name=angsdvcf
#SBATCH --ntasks=12
#SBATCH --nodes=3             
#SBATCH --time=40:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.angsdvcf.%j

module load parallel

angsting () { 

cd /xdisk/mcnew/dannyjackson/finches/vcf_likelihoods/"$@"/post

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/"$@"_post_bams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doBcf 1 -doGlf 3 -nThreads 12

cd /xdisk/mcnew/dannyjackson/finches/vcf_likelihoods/"$@"/pre

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/"$@"_pre_bams.txt -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doBcf 1 -doGlf 3 -nThreads 12

}

export -f angsting 

parallel -j 3 angsting ::: cra for par

sbatch angsdvcfs.sh 
Submitted batch job 2192337



#!/bin/bash

#SBATCH --job-name=makevcf
#SBATCH --ntasks=3
#SBATCH --nodes=3             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.makevcf.%j

module load bcftools 
module load parallel 

makevcf () {

cd /xdisk/mcnew/dannyjackson/finches/vcf_likelihoods/"$@"/pre

bcftools convert -O z -o angsdput.vcf.gz angsdput.bcf

cd /xdisk/mcnew/dannyjackson/finches/vcf_likelihoods/"$@"/post

bcftools convert -O z -o angsdput.vcf.gz angsdput.bcf


}

export -f makevcf 

parallel -j 3 makevcf ::: cra for par

sbatch makevcf.sh 
Submitted batch job 2192768
