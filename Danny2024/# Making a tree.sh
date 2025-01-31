# Making a tree

#!/bin/bash

#SBATCH --job-name=dobcf
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=32gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.dobcf.%j

module load bcftools 
cd /xdisk/mcnew/dannyjackson/finches/angsty/analyses/raxml

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt -out allsamples.bcf -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -doBcf 1

bcftools convert -O z -o allsamples.vcf.gz allsamples.bcf

sbatch dobcf.sh 
Submitted batch job 2163092




zcat allsamples.vcf.gz | grep 'CHROM' | head -1

gunzip allsamples.vcf.gz

sed -i 's\/xdisk\/mcnew\/dannyjackson\/finches\/bias\_testing\/batchnaive\/indelrealignment\///g' allsamples.vcf

sed 's/xxx/yyy/g'
sed -i 's/\.realigned\.bam//g' allsamples.vcf




