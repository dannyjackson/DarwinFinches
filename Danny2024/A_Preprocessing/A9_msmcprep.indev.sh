# I don't think I need to do any of this since the cardinal msmc prep processed the same reference genome

# msmc
cp -r /xdisk/mcnew/dannyjackson/cardinals/analyses/msmc/msmc2_scripts .
cd msmc2_scripts/

cd ~/programs
git clone https://github.com/stschiff/msmc2
git clone https://github.com/stschiff/msmc-tools
git clone https://github.com/jessicarick/msmc2_scripts

cd ~/programs/msmc2_scripts/
# edit msmc_params.sh
pip3 install --user whatshap

# following Jessi Rick's tutorial, found here: https://github.com/jessicarick/msmc2_scripts?tab=readme-ov-file
# Step 0: Create mappability mask -- already completed for this reference genome in my cardinalis analyses
# Step 1: generate vcf and mask files for each individual and each chromosome


cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/msmc

# add these lines to msmc_1_call.sh and comment out previous versions:
BAMFILE=$1
IND=$(basename "$BAMFILE" .realigned.bam)  # Extract filename without path and extension


scriptdir=$(dirname "$0")

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt 
~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/


# A2.4 individual mask vcf 
# test
sbatch --account=mcnew \
--job-name=ind_mask_vcf \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.ind_mask_vcf.%j \
--nodes=1 \
--ntasks-per-node=8 \
--time=48:00:00 \
~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ~/programs/DarwinFinches/params_preprocessing.sh -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ -i JP4481_all
# 3925911

# run in a slurm array
for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt `;
	do echo $i
      IND=$(basename "$i" .realigned.bam)  # Extract filename without path and extension

    sbatch --account=mcnew \
    --job-name=mask_vcf_${IND} \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.mask_vcf_${IND}.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ~/programs/DarwinFinches/params_preprocessing.sh -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ -i $IND
done

# do for duplicated bams

# run in a slurm array
for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/duplicatedsamples.txt `;
	do echo $IND
    sbatch --account=mcnew \
    --job-name=mask_vcf_${IND} \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.mask_vcf_${IND}.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.4_individual_mask_vcf.sh -p ~/programs/DarwinFinches/params_preprocessing.sh -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ -i $IND
done
12475825..12475848


# phasing

for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt `; 
  do echo $i 
    IND=$(basename "$i" .realigned.bam)  # Extract filename without path and extension

  sbatch --account=mcnew \
  --job-name=phase_${IND} \
  --partition=standard \
  --mail-type=ALL \
  --output=slurm_output/output.phase_${IND}.%j \
  --nodes=1 \
  --ntasks-per-node=8 \
  --time=48:00:00 \
  ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_phasing.sh -p ~/programs/DarwinFinches/params_preprocessing.sh -b /xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/ -i $IND
done

/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf/JP4481_all.NC_044571.1.vcf

3926689..3926739