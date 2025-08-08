# Indel Realignment Pipeline – GATK3

This document contains annotated SLURM job scripts for performing **Indel Realignment** using **GATK3** inside an Apptainer container. The workflow has two main steps:

1. **RealignerTargetCreator** – finds intervals likely to benefit from realignment.
2. **IndelRealigner** – realigns reads around those intervals.

**Goal:** For each *duplicate-marked, coordinate-sorted, overlap-clipped* BAM (`*.all.sorted.marked.clipped.bam`), realign around indels.

**Tool:** [`gatk3_3.7-0 GenomeAnalysisTK`] (invoked via apptainer exec)
---

## Step 1 – RealignerTargetCreator

```bash
## Indel realignment

#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools

# Loop over each sample name from the list, then index BAM file for GATK, then run RealignerTargetCreator to find indel intervals

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

---

## Step 2 – IndelRealigner

```bash
#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j

# Loop over sample names and realign reads to identified intervals
while read -r finch;
do 
  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
  -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals \
  -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

---

## Running the Jobs

```bash
sbatch indelmap.sh
# Submits the RealignerTargetCreator job

sbatch indelrealignment.sh
# Submits the IndelRealigner job (be sure the first is donw before running this one)
```

---

## Alternate Version – `.all` Suffix Samples

This is the same as above, but operates on BAM files with `.all.sorted.marked.clipped.bam` suffix and outputs `*_all` filenames. This is necessary because some of our samples were sequenced in multiple runs which were combined into "_all" files.

**RealignerTargetCreator – All Samples**

```bash
#!/bin/bash

#SBATCH --job-name=indelmap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelmap.%j

module load samtools

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch"_all.intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt 
```

**IndelRealigner – All Samples**

```bash
#!/bin/bash

#SBATCH --job-name=indelrealignment
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.indelrealignment.%j

while read -r finch;
do 
  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar -T IndelRealigner \
  -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
  --consensusDeterminationModel USE_READS \
  -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".all.sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch"_.intervals \
  -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch"_all.realigned.bam

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt 
```
---

## Running the Jobs

```bash
sbatch ~/programs/slurmscripts/indelrealignment.all.sh

sbatch ~/programs/slurmscripts/indelmap.all.sh 
```