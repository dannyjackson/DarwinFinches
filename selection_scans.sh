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

pixy --stats pi fst dxy --vcf ${vcf} --populations ${pixypop} --window_size ${windowsize} --output_folder ${outDir}/pixy --bypass_invariant_check yes

fi

