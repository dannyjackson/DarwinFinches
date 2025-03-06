# msmc

cd ~/programs
git clone https://github.com/stschiff/msmc2
git clone https://github.com/stschiff/msmc-tools
git clone https://github.com/jessicarick/msmc2_scripts

cd ~/programs/msmc2_scripts/
# edit msmc_params.sh
pip3 install --user whatshap

# following Jessi Rick's tutorial, found here: https://github.com/jessicarick/msmc2_scripts?tab=readme-ov-file
# Step 0: Create mappability mask

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/msmc


sbatch --account=mcnew \
--job-name=run_snpable2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.run_snpable2.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
/xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/msmc/msmc2_scripts/run_snpable2.sh

# Submitted batch job 11143339

awk '{print $1}' /xdisk/mcnew/finches/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa.fai > /xdisk/mcnew/dannyjackson/cardinals_dfinch/analyses/msmc/msmc2_scripts/SCAFFOLDS.txt

# then edit the makeMappabilityMask.py script (lines 26 & 30) for your genome, and run the script
# on line 30, use curly braces {} to indicate where in the name the scaffold name should go

# /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

splitfa $GENOME $k | split -l 20000000
cat snpable/x* >> GCF_963921805__split.150

gzip -dc xx??.sam.gz | gen_raw_mask.pl > rawMask_35.fa  

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${GENOME} ${prefix}_split.${k} > ${prefix}_split.${k}.sai
bwa samse -f ${prefix}_split.${k}.sam $GENOME ${prefix}_split.${k}.sai ${prefix}_split.${k}


sbatch --account=mcnew \
--job-name=snpable_2 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.snpable_2.%j \
--nodes=1 \
--ntasks-per-node=16 \
--time=48:00:00 \
snpable_2.sh

Submitted batch job 3526455


#!/bin/bash
module load bwa/0.7.17
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module spider samtools/1.19.2

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/msmc

bwa aln -t 8 -R 1000000 -O 3 -E 3 /xdisk/mcnew/dannyjackson/sulidae/old/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna GCF_963921805__split.150 > GCF_963921805__split.150.sai

bwa samse -f GCF_963921805__split.150.sam /xdisk/mcnew/dannyjackson/sulidae/old/angsd/refgenome/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna GCF_963921805__split.150.sai GCF_963921805__split.150 

sbatch snpable_2.sh 
Submitted batch job 3579511