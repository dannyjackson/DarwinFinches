#!/bin/bash


# shell script to pull significant genes from an analysis of a vcf containing two populations: pop1 is of interest and pop2 is the reference

# for example, when looking for genes under selection in an urban population, pop1 would be the urban samples and pop2 would be the rural samples

# analyzes sliding windows of FST, nucleotide diversity (pi) and Tajima's D.

#
#

if [ $# -lt 1 ]
  then
    echo "Analyzes a vcf containing sequences from two populations using fst, nucleotide diversity, and Tajima's D,.

    [-v] Path to full vcf file 
    [-n] Output project name (don't include extension)
    [-o] Output directory for files
    [-p] Path to pop1 vcf (reference population)
    [-q] Path to pop2 vcf (population of interest)
    [-r] Path to pixy population list file
    [-g] Path to gff file"


  else
    while getopts v:n:o:p:q:r:g: option
    do
    case "${option}"
    in
    v) vcf=${OPTARG};;
    n) name=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pop1=${OPTARG};;
    q) pop2=${OPTARG};;
    r) pixypop=${OPTARG};;
    g) gff=${OPTARG};;
    esac
    done


module load vcftools
module load R



##  fst, pi, and dxy
mkdir ${outDir}/pixy

source activate pixy

module load samtools

bgzip ${vcf}
tabix ${vcf}.gz

pixy --stats fst pi dxy --vcf ${vcf} --populations ${pixypop} --window_size 1000 --output_folder ${outDir}/pixy --bypass_invariant_check yes

head -1 ${outDir}/pixy/pixy_dxy.txt > ${outDir}/pixy/${name}.chroms.pixy_dxy.txt
grep 'NC' ${outDir}/pixy/pixy_dxy.txt >> ${outDir}/pixy/${name}.chroms.pixy_dxy.txt


sed -i 's/NC_//g' ${outDir}/pixy/${name}.chroms.pixy_dxy.txt

awk '{sub(/\./,"",$1)}1' ${outDir}/pixy/${name}.chroms.pixy_dxy.txt | column -t > ${outDir}/pixy/${name}.chroms.pixy_dxy.formanhattan.txt

Rscript ~/programs/DarwinFinches/pixy.r ${outDir}/pixy/ ${name}

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/referencepop/${name}.tD_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/referencepop/${name}.tD_sig_genes.csv

## Tajima's D

# population of interest

mkdir ${outDir}/tajimasd 
mkdir ${outDir}/tajimasd/referencepop ${outDir}/tajimasd/interestpop

vcftools --vcf ${pop2} --TajimaD 20000 --out ${outDir}/tajimasd/interestpop/${name} 


head -1 ${outDir}/tajimasd/interestpop/${name}.Tajima.D > ${outDir}/tajimasd/interestpop/${name}.chroms.Tajima.D
grep 'NC' ${outDir}/tajimasd/interestpop/${name}.Tajima.D >> ${outDir}/tajimasd/interestpop/${name}.chroms.Tajima.D


sed -i 's/NC_//g' ${outDir}/tajimasd/interestpop/${name}.chroms.Tajima.D

awk '{sub(/\./,"",$1)}1' ${outDir}/tajimasd/interestpop/${name}.chroms.Tajima.D | column -t > ${outDir}/tajimasd/interestpop/${name}.chroms.Tajima.D.formanhattan

Rscript ~/programs/DarwinFinches/tajimasD.r ${outDir}/tajimasd/interestpop ${name}

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/interestpop/${name}.tD_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/interestpop/${name}.tD_sig_genes.csv

# reference population

vcftools --vcf ${pop1} --TajimaD 20000 --out ${outDir}/tajimasd/referencepop/${name} 


head -1 ${outDir}/tajimasd/referencepop/${name}.Tajima.D > ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D
grep 'NC' ${outDir}/tajimasd/referencepop/${name}.Tajima.D >> ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D


sed -i 's/NC_//g' ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D

awk '{sub(/\./,"",$1)}1' ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D | column -t > ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D.formanhattan

Rscript ~/programs/DarwinFinches/tajimasD.r ${outDir}/tajimasd/referencepop ${name}

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/referencepop/${name}.tD_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/referencepop/${name}.tD_sig_genes.csv


fi

