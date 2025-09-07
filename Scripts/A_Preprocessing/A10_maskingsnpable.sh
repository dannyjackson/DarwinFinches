#!/bin/bash
# Genotype Call Preprocessing Pipeline for Darwin Finches
# 1. Generate VCF and Mask Files for Each Individual
# 2. Phasing
# 3. Sort and Index All VCF Files

############################################
# 0. Clone Required Repositories & Install Dependencies
############################################

cd ~/programs
git clone https://github.com/stschiff/msmc2
git clone https://github.com/stschiff/msmc-tools
git clone https://github.com/jessicarick/msmc2_scripts

cd ~/programs/msmc2_scripts/
pip3 install --user whatshap

# Reference: https://github.com/jessicarick/msmc2_scripts?tab=readme-ov-file

############################################
# 1. Generate VCF and Mask Files for Each Individual
############################################

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/msmc

# Run individual mask VCF generation for all samples
while read BAMFILE; do
    IND=$(basename "$BAMFILE" .realigned.bam)
    sbatch --account=mcnew \
        --job-name=mask_vcf_${IND} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.mask_vcf_${IND}.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=48:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh \
        -p ~/programs/DarwinFinches/params_preprocessing.sh \
        -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ \
        -i $IND
done < /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt

# Repeat for duplicated samples
while read IND; do
    sbatch --account=mcnew \
        --job-name=mask_vcf_${IND} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.mask_vcf_${IND}.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=48:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh \
        -p ~/programs/DarwinFinches/params_preprocessing.sh \
        -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ \
        -i $IND
done < /xdisk/mcnew/finches/dannyjackson/finches/referencelists/duplicatedsamples.txt

############################################
# 2. Phasing
############################################

while read BAMFILE; do
    IND=$(basename "$BAMFILE" .realigned.bam)
    sbatch --account=mcnew \
        --job-name=phase_${IND} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.phase_${IND}.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=48:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_phasing.sh \
        -p ~/programs/DarwinFinches/params_preprocessing.sh \
        -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ \
        -i $IND
done < /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt

# Repeat for duplicated samples
while read IND; do
    sbatch --account=mcnew \
        --job-name=phase_${IND} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.phase_${IND}.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=48:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_phasing.sh \
        -p ~/programs/DarwinFinches/params_preprocessing.sh \
        -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ \
        -i $IND
done < /xdisk/mcnew/finches/dannyjackson/finches/referencelists/duplicatedsamples.txt

############################################
# 3. Sort and Index All VCF Files
############################################

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2
find . -name "*.vcf.gz" > list.txt

module load parallel
module load bcftools

cat list.txt | parallel -j 12 '
    file={}
    sorted_file=${file/input_file/sorted_input_file}
    bcftools sort "$file" -Oz -o "$sorted_file"
    bcftools index -t "$sorted_file"
'

# Submit sorting/indexing job via SLURM
sbatch --account=mcnew \
    --job-name=sortindex \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.sortindex.%j \
    --nodes=1 \
    --ntasks-per-node=12 \
    --time=48:00:00 \
    /xdisk/mcnew/finches/dannyjackson/finches/analyses/msmc/sort_index_vcf2.sh
