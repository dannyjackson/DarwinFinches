#!/bin/bash


# shell script to pull significant genes from an analysis of a vcf containing two populations: pop1 is of interest and pop2 is the reference

# for example, when looking for genes under selection in an urban population, pop1 would be the urban samples and pop2 would be the rural samples

# uses  

#
#

if [ $# -lt 1 ]
  then
    echo "Analyzes a vcf containing sequences from two populations using fst, bayescan, nucleotide diversity, Tajima's D, and SweeD.

    [-v] Path to full vcf file 
    [-n] Output project name (don't include extension)
    [-o] Output directory for files
    [-p] Path to pop1 vcf (reference population)
    [-q] Path to pop2 vcf (population of interest)
    [-r] Path to pop1 sample list file (reference population)
    [-s] Path to pop2 sample list file (population of interest)
    [-g] Path to gff file"


  else
    while getopts v:n:o:p:q:r:s:g: option
    do
    case "${option}"
    in
    v) vcf=${OPTARG};;
    n) name=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pop1=${OPTARG};;
    q) pop2=${OPTARG};;
    r) p1file=${OPTARG};;
    s) p2file=${OPTARG};;
    g) gff=${OPTARG};;
    esac
    done


module load vcftools
module load R

mkdir ${outDir}/fst 

##  fst

vcftools --vcf ${vcf} --weir-fst-pop ${p1file} --weir-fst-pop ${p2file} --out ${outDir}/fst/${name} --fst-window-size 10000 


head -1 ${outDir}/fst/${name}.windowed.weir.fst > ${outDir}/fst/${name}.chroms.windowed.weir.fst
grep 'NC' ${outDir}/fst/${name}.windowed.weir.fst >> ${outDir}/fst/${name}.chroms.windowed.weir.fst

sed -i 's/NC_//g' ${outDir}/fst/${name}.chroms.windowed.weir.fst

Rscript ~/programs/DarwinFinches/fstscans.r ${outDir}/fst ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/fst/${name}.fst_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/fst/${name}.zfst_sig.csv

## bayescan

## nucleotide diversity
mkdir ${outDir}/nucleotidediversity 
mkdir ${outDir}/nucleotidediversity/referencepop ${outDir}/nucleotidediversity/interestpop

# reference population analysis
vcftools --vcf ${pop1} --window-pi 10000 --out ${outDir}/nucleotidediversity/referencepop/${name}

head -1 ${outDir}/nucleotidediversity/referencepop/${name}.windowed.pi > ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi
grep 'NC' ${outDir}/nucleotidediversity/referencepop/${name}.windowed.pi >> ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi | column -t > ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi.formanhattan

Rscript ~/programs/DarwinFinches/nucleotidediversity.r ${outDir}/nucleotidediversity/referencepop ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/nucleotidediversity/referencepop/${name}.pi_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/nucleotidediversity/referencepop/${name}.pi_sig_genes.csv

# population of interest analysis
vcftools --vcf ${pop2} --window-pi 10000 --out ${outDir}/nucleotidediversity/interestpop/${name}

head -1 ${outDir}/nucleotidediversity/interestpop/${name}.windowed.pi > ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi
grep 'NC' ${outDir}/nucleotidediversity/interestpop/${name}.windowed.pi >> ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi | column -t > ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi.formanhattan

Rscript ~/programs/DarwinFinches/nucleotidediversity.r ${outDir}/nucleotidediversity/interestpop ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/nucleotidediversity/interestpop/${name}.pi_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/nucleotidediversity/interestpop/${name}.pi_sig_genes.csv

## Tajima's D
## SweeD

fi

