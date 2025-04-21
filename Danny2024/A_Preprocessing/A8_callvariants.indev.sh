#!/bin/bash
module load bcftools/1.19

ref="/xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
ID="finches"

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/

bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP --bam-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


sbatch --account=mcnew \
--job-name=callvariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.callvariants.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=100:00:00 \
callvariants.sh

# Submitted batch job 
# 100 hrs 3916091


#!/bin/bash

module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module spider samtools/1.19.2

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/

bcftools view -i 'QUAL>100' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.vcf > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_qualitysort.vcf

#filters by depth and removes indels
vcftools --vcf /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_qualitysort.vcf --min-meanDP 2 --max-meanDP 8 --remove-indels --recode --out /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_filtered

plink --vcf /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_filtered_mind2


sbatch --account=mcnew \
--job-name=filtervariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.filtervariants.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=100:00:00 \
filtervariants.sh

# Submitted batch job 12491252