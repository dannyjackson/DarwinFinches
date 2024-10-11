#!/bin/bash

#SBATCH --job-name=realign_fixmatess
#SBATCH --ntasks=94
#SBATCH --nodes=1             
#SBATCH --time=72:01:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mcnew@arizona.edu
#SBATCH --output=output.realign_fixmates.%j
#SBATCH --constraint=hi_mem

# update from Cornell script
# now need to move to gatk 4, and follow best practices. 
# no longer necessary to do Indel realigner or target creator 
module load picard
module load samtools 
module load gatk
module list

cd /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs

# List of all samples 
#finch=(JP4481 JP5410 JP9655 lamich_PARV1 lamich_PARV2 lamich_PL15 lamich_PL16 lamich_PL4 lamich_PL7 lamich_PL9 RHC097 RHC507 SM031 SM032 SM040 SM059 SM079 SM1067 SM1083 SM1156 SM1157 SM1200 SM1204 SM1231 SM1237 SM1240 SM1266 SM1270 SM1271 SM1272 SM1273 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338)

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507 JP4481_round1 JP9655_round1 RHC507_round1 SM032_round1 SM059_round1 SM1067_round1 SM1200_round1 SM1240_round1 JP5410_round1 RHC097_round1 SM031_round1 SM040_round1 SM079_round1 SM1157_round1 SM1204_round1 SM1266_round1)

# combine Markduplicate metrics into file 
for u in "${finch[@]}"; do
cat ${u}.duplicate.metrics.txt >> mark.duplicate.metrics.txt
done

# sort 
for u in "${finch[@]}"; do
samtools index  ${u}.sorted.marked.bam
done

# Fix sample info for double-library samples so that Haplotype Caller works 
# Need to rename sample name in each of the lamich samples. Use Samtools to pull header
# then manually change the SM: field in the header. use. "Reheader" to add header back on. 
# sed -i 's/SM\:SRR1607505/SM\:lamich_\PARV1/g' header etc
# Resources: https://gatk.broadinstitute.org/hc/en-us/articles/360035889471-How-should-I-pre-process-data-from-multiplexed-sequencing-and-multi-library-designs-
# https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/frequently-asked-questions/3060-How-should-I-preprocess-data-from-multiplexed-sequencing-and-multilibrary-designs

finch=(lamich_PARV1 lamich_PARV2 lamich_PL15 lamich_PL16 lamich_PL4 lamich_PL7 lamich_PL9)

for u in "${finch[@]}"; do
samtools view -H ${u}.sorted.marked.bam > header.${u}.sam
done

for u in "${finch[@]}"; do
samtools reheader header.${u}.sam ${u}.sorted.marked.bam > ${u}.sorted.marked.rehead.bam
done

for u in "${finch[@]}"; do
samtools index ${u}.sorted.marked.rehead.bam
done




echo "Job ended on `date`"
