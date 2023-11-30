#!/bin/bash


# shell script to pull significant genes from an analysis of a vcf containing two populations: pop1 is of interest and pop2 is the reference

# for example, when looking for genes under selection in an urban population, pop1 would be the urban samples and pop2 would be the rural samples

# uses  

#
#

if [ $# -lt 1 ]
  then
    echo "Analyzes a vcf containing sequences from two populations using fst, bayescan, nucleotide diversity, Tajima's D, and SweeD.

    [-i] Path to vcf file directory
    [-v] Vcf file name (don't include extension)
    [-o] Output directory for files
    [-p] Path to pop1 file (population of interest)
    [-q] Path to pop2 file (reference population)"


  else
    while getopts i:o:c:r: option
    do
    case "${option}"
    in
    i) inDir=${OPTARG};;
    v) vcf=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pop1=${OPTARG};;
    q) pop2=${OPTARG};;

    esac
    done


module load vcftools
module load R


vcftools --vcf ${inDir}/${vcf}.vcf --window-pi 10000 --out ${outDir}/${vcf}

head -1 ${outDir}/${vcf}.windowed.pi > ${outDir}/${vcf}.chroms.windowed.pi
grep 'NC' ${outDir}/${vcf}.windowed.pi >> ${outDir}/${vcf}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/${vcf}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/${vcf}.chroms.windowed.pi | column -t > ${outDir}/${vcf}.chroms.windowed.pi.formanhattan



Rscript nucleotidediversity.R ${outDir} ${vcf}

fi