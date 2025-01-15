#!/bin/bash

#SBATCH --job-name=trim_fastqs
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=24:01:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.trim_fastqs.%j

# Script to remove adaptors 
cd /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/
export PATH=/groups/mcnew/software/adapterremoval-2.3.1/build/:$PATH

echo "Job started on `date`"

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)

for u in "${finch[@]}"; do
AdapterRemoval --file1 ${u}_1.fastq.gz --file2 ${u}_2.fastq.gz --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 8  --basename ${u}
done


finch=(JP4481 JP9655 RHC507 SM032 SM059 SM1067 SM1200 SM1240 JP4481 JP9655 RHC507 SM032 SM059 SM1067 SM1200 SM1240 JP5410 RHC097 SM031 SM040 SM079 SM1157 SM1204 SM1266 JP5410 RHC097 SM031 SM040 SM079 SM1157 SM1204 SM1266)

for u in "${finch[@]}"; do
AdapterRemoval --file1 ${u}_1.round1.fq.gz --file2 ${u}_2.round1.fq.gz  --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 8  --basename ${u}_round1
done

echo "Job ended on `date`"
