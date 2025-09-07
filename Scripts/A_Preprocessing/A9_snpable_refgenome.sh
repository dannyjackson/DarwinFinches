#!/bin/bash
# Genotype Call Preprocessing Pipeline for Darwin Finches
# Mapping callable regions of the genome with SNPable

source ~/programs/DarwinFinches/param_files/params_base.sh

snpable_script_path=~/programs/seqbility-20091110 # directory with snpable scripts
PATH=$PATH:$snpable_script_path

date

scriptdir=${OUTDIR}/analyses/msmc/msmc2_scripts
source ${scriptdir}/msmc_params.sh

mkdir ${OUTDIR}/snpable
cd ${OUTDIR}/snpable

echo "Starting extraction of overlapping ${k}-mer subsequences"

splitfa $GENOME $k | split -l 20000000
cat x* >> ${prefix}_split.$k

# if it can't find splitfa, try adding seqbility to the path using 'PATH=$PATH:/project/WagnerLab/jrick/msmc_Sept2017/snpable/scripts'

echo "Aligning ${k}-mer reads to the genome with BWA, then converting to sam file"

# the genome needs to be indexed prior to this step-- if it has not already been indexed, run:
if [ -f "${GENOME}.bwt" ]; then
        echo "$GENOME already indexed"
else
        echo "indexing $GENOME"
        bwa index $GENOME
fi

echo "aligning reads to genome with BWA and converting to sam"
bwa aln -t 8 -R 1000000 -O 3 -E 3 ${GENOME} ${prefix}_split.${k} > ${prefix}_split.${k}.sai
bwa samse -f ${prefix}_split.${k}.sam $GENOME ${prefix}_split.${k}.sai ${prefix}_split.${k}

echo "reads aligned, starting to generate rawMask"
gen_raw_mask.pl ${prefix}_split.${k}.sam > ${prefix}_rawMask.${k}.fa

echo "raw mask created as ${prefix}_rawMask.35.fa, now generating final mask with stringency r=50%"
gen_mask -l ${k} -r 0.5 ${prefix}_rawMask.${k}.fa > ${prefix}_mask.${k}.50.fa

echo "all done! final mask saved as ${prefix}_mask.${k}.50.fa"

date
