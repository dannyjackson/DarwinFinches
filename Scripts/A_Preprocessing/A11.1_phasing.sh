#!/bin/bash

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
			whatshap stats --tsv=${OUTDIR}/stats/${IND}.${s}.${prefix}.minDP10.${phasing}.stats.tsv ${OUTDIR}/vcf3/${IND}.${s}.${phasing}.samtools.vcf.gz ;
    fi;
done

echo "finished with individual $IND"
date