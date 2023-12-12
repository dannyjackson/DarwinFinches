#!/bin/bash

#SBATCH --job-name=selection_scans
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.selection.%j


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

pixy --stats pi fst dxy --vcf ${vcf}.gz --populations ${pixypop} --window_size ${windowsize} --output_folder ${outDir}/pixy --bypass_invariant_check yes


mkdir ${outDir}/pixy/pi ${outDir}/pixy/fst ${outDir}/pixy/dxy 

mv ${outDir}/pixy/pixy_pi.txt ${outDir}/pixy/pi
mv ${outDir}/pixy/pixy_fst.txt ${outDir}/pixy/fst
mv ${outDir}/pixy/pixy_dxy.txt ${outDir}/pixy/dxy

## preparing files for r script

## fst

head -1 ${outDir}/pixy/fst/pixy_fst.txt > ${outDir}/pixy/fst/${name}.chroms.pixy_fst.txt
grep 'NC' ${outDir}/pixy/fst/pixy_fst.txt >> ${outDir}/pixy/fst/${name}.chroms.pixy_fst.txt


sed -i 's/NC_//g' ${outDir}/pixy/fst/${name}.chroms.pixy_fst.txt

awk '{sub(/\./,"",$1)}1' ${outDir}/pixy/fst/${name}.chroms.pixy_fst.txt | column -t > ${outDir}/pixy/fst/${name}.chroms.pixy_fst.formanhattan.txt


## pi

head -1 ${outDir}/pixy/pi/pixy_pi.txt > ${outDir}/pixy/pi/${name}.chroms.pixy_pi.txt
grep 'NC' ${outDir}/pixy/pi/pixy_pi.txt >> ${outDir}/pixy/pi/${name}.chroms.pixy_pi.txt


sed -i 's/NC_//g' ${outDir}/pixy/pi/${name}.chroms.pixy_pi.txt

awk '{sub(/\./,"",$1)}1' ${outDir}/pixy/pi/${name}.chroms.pixy_pi.txt | column -t > ${outDir}/pixy/pi/${name}.chroms.pixy_pi.formanhattan.txt


## dxy
head -1 ${outDir}/pixy/dxy/pixy_dxy.txt > ${outDir}/pixy/dxy/${name}.chroms.pixy_dxy.txt
grep 'NC' ${outDir}/pixy/dxy/pixy_dxy.txt >> ${outDir}/pixy/dxy/${name}.chroms.pixy_dxy.txt


sed -i 's/NC_//g' ${outDir}/pixy/dxy/${name}.chroms.pixy_dxy.txt

awk '{sub(/\./,"",$1)}1' ${outDir}/pixy/dxy/${name}.chroms.pixy_dxy.txt | column -t > ${outDir}/pixy/dxy/${name}.chroms.pixy_dxy.formanhattan.txt

Rscript ~/programs/DarwinFinches/pixy.r ${outDir}/pixy/ ${name}



## post R script 
## fst 
awk -F"[,\t]" 'NR==FNR{a["NC_0"$4]=$0; b=($5-10000); c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/pixy/fst/${name}.zfst_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/pixy/fst/${name}.zfst_sig_genes.csv

## dxy
awk -F"[,\t]" 'NR==FNR{a["NC_0"$4]=$0; b=($5-10000); c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/pixy/dxy/${name}.zdxy_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/pixy/fst/${name}.dxy_sig_genes.csv


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

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/interestpop/${name}.tD_sig.tsv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/interestpop/${name}.tD_sig_genes.tsv

# reference population

vcftools --vcf ${pop1} --TajimaD 20000 --out ${outDir}/tajimasd/referencepop/${name} 


head -1 ${outDir}/tajimasd/referencepop/${name}.Tajima.D > ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D
grep 'NC' ${outDir}/tajimasd/referencepop/${name}.Tajima.D >> ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D


sed -i 's/NC_//g' ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D

awk '{sub(/\./,"",$1)}1' ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D | column -t > ${outDir}/tajimasd/referencepop/${name}.chroms.Tajima.D.formanhattan

Rscript ~/programs/DarwinFinches/tajimasD.r ${outDir}/tajimasd/referencepop ${name}

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/referencepop/${name}.tD_sig.tsv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/referencepop/${name}.tD_sig_genes.tsv

fi

