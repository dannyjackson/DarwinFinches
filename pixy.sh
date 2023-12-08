#!/bin/bash

#SBATCH --job-name=pixy_cra
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pixy_cra.%j


if [ $# -lt 1 ]
  then
    echo "Analyzes a vcf containing sequences from two populations using fst, pi, and dxy.

    [-v] Path to full vcf file 
    [-o] Output directory for files
    [-r] Path to pixy pop file
    [-w] Window size"
 


  else
    while getopts v:o:r:w: option
    do
    case "${option}"
    in
    v) vcf=${OPTARG};;
    o) outDir=${OPTARG};;
    r) pixypop=${OPTARG};;
    w) windowsize=${OPTARG};;
    esac
    done

conda activate pixy

module load samtools

pixy --stats pi fst dxy --vcf ${vcf} --populations ${pixypop} --window_size ${windowsize} --output_folder ${outDir} --bypass_invariant_check yes

fi