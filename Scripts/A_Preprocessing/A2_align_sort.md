# Align and Sort Finch WGS Reads

This workflow uses **Bowtie2** for alignment, **Samtools** for BAM conversion & sorting, and **SLURM** for job scheduling.
It assumes you have **paired-end fastq.gz files** that have been trimmed, and a Bowtie2 index already built from the reference genome.

---

## 1. Identify Sample Names

```bash
# List all files containing "_2P" (paired-end read 2 after trimming)
# Extract sample name (everything before the first "_")
# Sort, remove duplicates, and count
ls *_2P* | awk 'BEGIN {FS = "_"} {print $1}' | sort | uniq -c
```

**Purpose:**

* Ensures you have a unique list of sample names before looping through them in alignment.

---

## 2. SLURM Job Script for Alignment

```bash
#!/bin/bash
#SBATCH --job-name=align_sequences
#SBATCH --ntasks=24
#SBATCH --nodes=1
#SBATCH --time=240:00:00
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

# Load required modules
module load bowtie2
module load picard
module load samtools
module list

# Move to directory containing trimmed FASTQ files
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

echo "Job started on $(date)"

# Reference genome preparation (run once)
# wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Camarhynchus_parvulus/latest_assembly_versions/GCF_901933205.1_STF_HiC/GCF_901933205.1_STF_HiC_genomic.fna.gz
# gzip -d GCF_901933205.1_STF_HiC_genomic.fna.gz
# mv GCF_901933205.1_STF_HiC_genomic.fna GCF_901933205.fa
# bowtie2-build -f GCF_901933205.fa GCF_901933205.fa
# picard CreateSequenceDictionary R=GCF_901933205.fa O=GCF_901933205.fa.dict
# samtools faidx GCF_901933205.fa

# Array of finch sample names
finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)

echo "Aligning FASTQ files"

# Step 1: Align reads with Bowtie2
for u in "${finch[@]}"; do
    echo "Aligning sample: ${u}"
    bowtie2 -p 24 --phred33 --very-sensitive-local \
        -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
        -I 149 -X 900 \
        --rg-id ${u} --rg SM:${u} \
        -1 ${u}_trimmed2_1P.fq.gz \
        -2 ${u}_trimmed2_2P.fq.gz \
        -U ${u}_trimmed2_1U.fq.gz,${u}_trimmed2_2U.fq.gz \
        -S /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam \
        2>> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/alignment.output.txt
done

# Step 2: Convert SAM to BAM and sort
echo "Converting SAMs to BAMs and sorting"
for u in "${finch[@]}"; do
    samtools view -S -b \
        /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam \
        > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam

    samtools sort \
        /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam \
        -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done

# Step 3: Compute alignment statistics
echo "Computing stats on sorted BAMs"
for u in "${finch[@]}"; do
    echo ${u}
    samtools flagstat \
        /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt

echo "Computing stats on unsorted BAMs"
for u in "${finch[@]}"; do
    echo ${u}
    samtools flagstat \
        /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/alignmentstats.txt

echo "Job ended on $(date)"
```

---

## 3. Extracting Sample Names from BAM Files

```bash
# List sorted marked BAM files and extract sample names (before first '.')
ls /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles \
    | awk 'BEGIN {FS = "."} {print $1}' \
    > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

# For files matching '*all*bam', extract sample names
ls *all*bam \
    | awk 'BEGIN {FS = "."} {print $1}' \
    > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt
```

**Purpose:**

* Maintains reference lists of sample names for downstream processing.