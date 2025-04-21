# Phase data

# multisample phasing
/xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.vcf

#!/bin/sh 
module load python

scriptdir=/xdisk/mcnew/finches/dannyjackson/finches/analyses/msmc/msmc2_scripts
source ${scriptdir}/msmc_params.sh

echo "working with individual $IND"

for s in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/SCAFFOLDS.txt`;
    do echo "working with scaffold $s";
    if [ -f ${OUTDIR}/vcf2/${IND}.${s}.${phasing}.samtools.vcf.gz ]; then
            echo "phased VCF already exists; moving onto next scaffold";
    else
            echo "phased VCF does not exist; phasing VCF for scaffold $s"
            sed -i 's/^ //g' ${OUTDIR}/vcf2/${IND}.${s}.samtools.vcf

            whatshap phase --reference ${GENOME} --ignore-read-groups -o ${OUTDIR}/vcf3/${IND}.${s}.${phasing}.samtools.vcf.gz ${OUTDIR}/vcf2/${IND}.${s}.samtools.vcf ${BAMFILE}
            whatshap stats --tsv=${OUTDIR}/stats/${IND}.${s}.${prefix}.minDP10.${phasing}.stats.tsv ${OUTDIR}/vcf/${IND}.${s}.${phasing}.samtools.vcf.gz ;
    fi;
done

echo "finished with individual $IND"
date


for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsamplebams.txt`;

for i in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/duplicatedsamples.txt`;

# for i in `cat msmc2_scripts/test_samplebam.txt`;
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
	phasing.sh $i
done

3925855..3925905



# 11152196 - 11152229
for job in {3925514..3925564}; do
    scancel $job
done
