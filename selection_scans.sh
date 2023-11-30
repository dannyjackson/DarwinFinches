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
    [-n] Output project name (don't include extension)
    [-o] Output directory for files
    [-p] Pop1 file name (reference population)
    [-q] Pop2 file name (population of interest)"


  else
    while getopts n:o:p:q: option
    do
    case "${option}"
    in
    n) name=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pop1=${OPTARG};;
    q) pop2=${OPTARG};;

    esac
    done


module load vcftools
module load R

mkdir ${outDir}/referencepop ${outDir}/interestpop

# reference population analysis
vcftools --vcf ${pop1} --window-pi 10000 --out ${outDir}/referencepop/${name}

head -1 ${outDir}/referencepop/${name}.windowed.pi > ${outDir}/referencepop/${name}.chroms.windowed.pi
grep 'NC' ${outDir}/referencepop/${name}.windowed.pi >> ${outDir}/referencepop/${name}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/referencepop/${name}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/referencepop/${name}.chroms.windowed.pi | column -t > ${outDir}/referencepop/${name}.chroms.windowed.pi.formanhattan



Rscript ~/programs/DarwinFinches/nucleotidediversity.r ${outDir}/referencepop ${name}


fi