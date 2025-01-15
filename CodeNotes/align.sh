#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=94
#SBATCH --nodes=1             
#SBATCH --time=72:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

module load bowtie2
module load picard
module load samtools
module list
 
cd /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs

echo "Job started on `date`"
# Download genome 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Camarhynchus_parvulus/latest_assembly_versions/GCF_901933205.1_STF_HiC/GCF_901933205.1_STF_HiC_genomic.fna.gz
#gzip -d GCF_901933205.1_STF_HiC_genomic.fna.gz 
#mv GCF_901933205.1_STF_HiC_genomic.fna GCF_901933205.fa
#bowtie2-build -f GCF_901933205.fa GCF_901933205.fa
#picard CreateSequenceDictionary R= GCF_901933205.fa O= GCF_901933205.fa.dict

#samtools faidx GCF_901933205.fa

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507 JP4481_round1 JP9655_round1 RHC507_round1 SM032_round1 SM059_round1 SM1067_round1 SM1200_round1 SM1240_round1 JP4481_round1 JP9655_round1 RHC507_round1 SM032_round1 SM059_round1 SM1067_round1 SM1200_round1 SM1240_round1 JP5410_round1 RHC097_round1 SM031_round1 SM040_round1 SM079_round1 SM1157_round1 SM1204_round1 SM1266_round1 JP5410_round1 RHC097_round1 SM031_round1 SM040_round1 SM079_round1 SM1157_round1 SM1204_round1 SM1266)



for u in "${finch[@]}"; do
echo ${u}
bowtie2 -p 24 --phred33 --very-sensitive-local -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -I 149 -X 900 --rg-id ${u} --rg SM:${u} -1 ${u}.pair1.truncated -2 ${u}.pair2.truncated -U ${u}.collapsed.truncated, ${u}.singleton.truncated -S ${u}.sam 2>> alignment.output.txt
done

# convert to bams and sort
for u in "${finch[@]}"; do 
samtools view -S -b ${u}.sam > ${u}.bam
samtools sort ${u}.bam -o ${u}.sorted.bam 
done


# some alignment stats 
for u in "${finch[@]}"; do
echo ${u}
samtools flagstat ${u}.bam
done > alignmentstats.txt 

echo "Job ended on `date`"