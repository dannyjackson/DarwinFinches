#!/bin/bash
module load bcftools/1.19

ref="/xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
ID="finches"

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/

bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP -S /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


sbatch --account=mcnew \
--job-name=callvariants \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.callvariants.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=100:00:00 \
callvariants.sh

# Submitted batch job 3915974
# 100 hrs 3916084