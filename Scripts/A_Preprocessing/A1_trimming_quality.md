````markdown
# Trimming and Quality Assessment Pipeline

This document describes the steps used for trimming and quality assessment of whole-genome sequencing (WGS) reads from finch samples.  
It includes scripts for **FastQC**, **Trimmomatic**, **Poly-G trimming**, and post-trimming quality checks.

---

## 1. Raw Read Quality Assessment with FastQC

```bash
#!/bin/bash

#SBATCH --job-name=fastQC
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

```bash
module load fastqc/0.11.9
cd /xdisk/mcnew/finches/dannyjackson/finches/bias_testing/platform/raw_fastqcs
```

Load FastQC and set working directory for output.

```bash
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1.round1.fq.gz 
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_2.round1.fq.gz 
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1.fastq.gz
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_2.fastq.gz
```

Run **FastQC** on both forward (`*_1`) and reverse (`*_2`) reads, before trimming.

---

**Organizing Results:**

```bash
ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*fastqc.zip
mv /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*fastqc.zip .
scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/finches/bias_testing/platform/raw_fastqcs/ .
```

* List zipped FastQC reports
* Move to working directory
* Copy reports from HPC to local machine

**Unzip and Extract Reports:**

```bash
ls > filenames.txt
while read -r file; do unzip $file; done < filenames.txt
ls > filenames.txt
while read -r file; do cp "$file"/fastqc_report.html fastqc_reports/"$file"_fastqc_report.html; done < filenames.txt
```

* Create list of filenames
* Unzip each `.zip` file
* Copy HTML reports into `fastqc_reports/`

**Extract warnings and failures:**

```bash
grep 'warn' */fastqc_data.txt
grep 'fail' */fastqc_data.txt
```

Identify modules flagged as **WARN** or **FAIL** in FastQC output.

---

## 2. Trimming with Trimmomatic

```bash
#!/bin/bash

#SBATCH --job-name=trimming2
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

cd ~/programs/Trimmomatic
```

**Annotations:**

* Running **Trimmomatic** for adapter/quality trimming.

```bash
echo "Beginning trimming for "$ID >> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt

while read -r ID; do
  java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 12 \
    /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz \
    /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz \
    -baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/"$ID"_trimmed2.fq.gz \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90 >> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim2_log.txt
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt
```

* **PE mode**: Paired-end trimming
* Removes bases with quality < 20 from start/end, trims with sliding window, and enforces min length of 90 bp.

---

## 3. Poly-G Tail Trimming with Fastp

```bash
sed -i 's/\.round1//g' /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt
```

* Removes `.round1` from filenames in sample list.

```bash
#!/bin/bash

#SBATCH --job-name=polygtrimming
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

while read -r ID; do
  echo "Beginning polyg trimming for "$ID >> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 \
    -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_1P.fq.gz \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_2P.fq.gz \
    -O /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt
```

* Uses **fastp** to trim poly-G stretches common in Illumina NovaSeq reads.

---

## 4. Post-Trimming FastQC

**Adapter-trimmed reads:**

```bash
module load fastqc/0.11.9
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/
fastqc -t 12 *.fq.gz
```

**Poly-G-trimmed reads:**

```bash
module load fastqc/0.11.9
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/
fastqc -t 12 *.fastq.gz
```

**Extract warning/failure statistics:**

```bash
grep 'warn' */fastqc_data.txt > warnings.txt
grep 'fail' */fastqc_data.txt > failings.txt

# Generate summary statistics
grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_adaptertrimmed_fastas.txt
grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_adaptertrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' failings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_polygtrimmed_fastas.txt
awk 'BEGIN {FS = ">>"} {print $2}' warnings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_polygtrimmed_fastas.txt

grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_finaltrim_fastas.txt
grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_finaltrim_fastas.txt
```

* Summarizes which modules fail or warn after each trimming step.
* Useful for comparing read quality improvement over the pipeline.

---


```
