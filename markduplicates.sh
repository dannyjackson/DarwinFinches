#!/bin/bash

#SBATCH --job-name=markdups
#SBATCH --ntasks=94
#SBATCH --nodes=1             
#SBATCH --time=72:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.markdups.%j


module load picard
module load samtools 
module list

cd /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs

# List of all samples 
# finch=(JP4481 JP5410 JP9655 RHC097 RHC507 SM031 SM032 SM040 SM059 SM079 SM1067 SM1083 SM1156 SM1157 SM1200 SM1204 SM1231 SM1237 SM1240 SM1266 SM1270 SM1271 SM1272 SM1273 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607532 SRR1607533 SRR1607534 SRR1607535 SRR1607536 SRR1607537 SRR1607538 SRR1607539 SRR1607540 SRR1607541 SRR1607542 SRR1607504 SRR1607505 SRR1607506 SRR1607507)

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507 JP4481_round1 JP9655_round1 RHC507_round1 SM032_round1 SM059_round1 SM1067_round1 SM1200_round1 SM1240_round1 JP5410_round1 RHC097_round1 SM031_round1 SM040_round1 SM079_round1 SM1157_round1 SM1204_round1 SM1266_round1)
# Mark Duplicates: starts with sorted bam files 

for u in "${finch[@]}"; do
java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates INPUT= ${u}.sorted.bam OUTPUT= ${u}.sorted.marked.bam METRICS_FILE=${u}.duplicate.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
done 

# Mark Duplicates for samples that had multiple BAMS

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607532.sorted.bam I= SRR1607533.sorted.bam \
OUTPUT= lamich_PL15.sorted.marked.bam \
METRICS_FILE=lamich_PL15.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607534.sorted.bam I= SRR1607535.sorted.bam I= SRR1607536.sorted.bam \
OUTPUT= lamich_PL16.sorted.marked.bam \
METRICS_FILE=lamich_PL16.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607537.sorted.bam I= SRR1607538.sorted.bam \
OUTPUT= lamich_PL4.sorted.marked.bam \
METRICS_FILE=lamich_PL4.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607539.sorted.bam I= SRR1607540.sorted.bam \
OUTPUT= lamich_PL7.sorted.marked.bam \
METRICS_FILE=lamich_PL7.duplicate.metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607541.sorted.bam I= SRR1607542.sorted.bam \
OUTPUT= lamich_PL9.sorted.marked.bam \
METRICS_FILE=lamich_PL9.duplicate.metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607505.sorted.bam I= SRR1607504.sorted.bam \
OUTPUT= lamich_PARV1.sorted.marked.bam \
METRICS_FILE= lamich_PARV1.duplicate.metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SRR1607506.sorted.bam I= SRR1607507.sorted.bam \
OUTPUT= lamich_PARV2.sorted.marked.bam \
METRICS_FILE=lamich_PARV2.duplicate.metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000

# combining bams from multiple sequence runs
java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= JP4481.sorted.bam I= JP4481_round1.sorted.bam \
OUTPUT= JP4481.all.sorted.marked.bam \
METRICS_FILE=JP4481.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= RHC507_round1.sorted.bam I= RHC507.sorted.bam \
OUTPUT= RHC507.all.sorted.marked.bam \
METRICS_FILE=RHC507.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM059_round1.sorted.bam I= SM059.sorted.bam \
OUTPUT= SM059.all.sorted.marked.bam \
METRICS_FILE=SM059.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1200_round1.sorted.bam I= SM1200.sorted.bam \
OUTPUT= SM1200.all.sorted.marked.bam \
METRICS_FILE=SM1200.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= JP5410_round1.sorted.bam I= JP5410.sorted.bam \
OUTPUT= JP5410.all.sorted.marked.bam \
METRICS_FILE=JP5410.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM031_round1.sorted.bam I= SM031.sorted.bam \
OUTPUT= SM031.all.sorted.marked.bam \
METRICS_FILE=SM031.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM079_round1.sorted.bam I= SM079.sorted.bam \
OUTPUT= SM079.all.sorted.marked.bam \
METRICS_FILE=SM079.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1204_round1.sorted.bam I= SM1204.sorted.bam \
OUTPUT= SM1204.all.sorted.marked.bam \
METRICS_FILE=SM1204.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= JP9655_round1.sorted.bam I= JP9655.sorted.bam \
OUTPUT= JP9655.all.sorted.marked.bam \
METRICS_FILE=JP9655.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM032_round1.sorted.bam I= SM032.sorted.bam \
OUTPUT= SM032.all.sorted.marked.bam \
METRICS_FILE=SM032.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1067_round1.sorted.bam I= SM1067.sorted.bam \
OUTPUT= SM1067.all.sorted.marked.bam \
METRICS_FILE=SM1067.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1240_round1.sorted.bam I= SM1240.sorted.bam \
OUTPUT= SM1240.all.sorted.marked.bam \
METRICS_FILE=SM1240.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= RHC097_round1.sorted.bam I= RHC097.sorted.bam \
OUTPUT= RHC097.all.sorted.marked.bam \
METRICS_FILE=RHC097.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM040_round1.sorted.bam I= SM040.sorted.bam \
OUTPUT= SM040.all.sorted.marked.bam \
METRICS_FILE=SM040.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1157_round1.sorted.bam I= SM1157.sorted.bam \
OUTPUT= SM1157.all.sorted.marked.bam \
METRICS_FILE=SM1157.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= SM1266_round1.sorted.bam I= SM1266.sorted.bam \
OUTPUT= SM1266.all.sorted.marked.bam \
METRICS_FILE=SM1157.all.duplicate.metrics.txt
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &



echo "Job ended on `date`"
