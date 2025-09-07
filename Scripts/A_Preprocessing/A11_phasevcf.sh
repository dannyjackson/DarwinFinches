#!/bin/bash
# Genotype Call Preprocessing Pipeline for Darwin Finches
# Phase vcf files for each individual after masking 


for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt`;
	do echo $i
    IND=$(basename "$i" .realigned.bam)  # Extract filename without path and extension
    
	sbatch --account=mcnew \
	--job-name=phasing_${IND} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.phasing${IND}.%j \
	--nodes=1 \
	--ntasks-per-node=4 \
	--time=25:00:00 \
	A11.1_phasing.sh $i
done

for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/duplicatedsamples.txt`;

	do echo $i
    IND=$(basename "$i" .realigned.bam)  # Extract filename without path and extension
    
	sbatch --account=mcnew \
	--job-name=phasing_${IND} \
    --partition=standard \
	--mail-type=ALL \
	--output=slurm_output/output.phasing${IND}.%j \
	--nodes=1 \
	--ntasks-per-node=4 \
	--time=25:00:00 \
	A11.1_phasing.sh $i
done
