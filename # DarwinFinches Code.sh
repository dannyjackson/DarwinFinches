# DarwinFinches Code 

# NEXT STEPS:
1. Bayescan
2. Timesweeper
3. PSMC
4. dadi 

## fastqc

#!/bin/bash

#SBATCH --job-name=angsd_forpopgen2
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=60gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_forpopgen2.%j

module load fastqc/0.11.9

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/platform/raw_fastqcs


fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1.round1.fq.gz 
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_2.round1.fq.gz 


fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1.fastq.gz
fastqc -t 12 /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_2.fastq.gz

sbatch fastqc_rawfastas.sh 
Submitted batch job 1888141

ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*fastqc.zip

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/platform/raw_fastqcs
mv /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*fastqc.zip .

scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/finches/bias_testing/platform/raw_fastqcs/ .

ls > filenames.txt

while read -r file;
do 
 unzip $file
done < filenames.txt

ls > filenames.txt

while read -r file;
do 
 cp "$file"/fastqc_report.html fastqc_reports/"$file"_fastqc_report.html
done < filenames.txt

# output warn and fail flags from fastqc data
# unpack zip files and navigate to directory above all zip files
grep 'warn' */fastqc_data.txt
grep 'fail' */fastqc_data.txt

## Trimmomatic 


#!/bin/bash

#SBATCH --job-name=trimming2
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=trimming2.%j

cd ~/programs/Trimmomatic

echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt

while read -r ID;
do
  java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 12 \
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz  /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz  \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/"$ID"_trimmed2.fq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim2_log.txt


done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt


## polyg

sed -i 's/\.round1//g' /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

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

while read -r ID;
do

echo "Beginning polyg trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_1P.fq.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_2P.fq.gz -O /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz 

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt


sbatch polyg_trimming.sh 
Submitted batch job 9453675


## Post-trimming fastqc 
#!/bin/bash

#SBATCH --job-name=fastqc_trimmed
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=fastqc_trimmed.%j

module load fastqc/0.11.9

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/

fastqc -t 12 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/*.fq.gz

sbatch fastqc_adaptertrimmed.sh 
Submitted batch job 9456471

#!/bin/bash

#SBATCH --job-name=fastqc_trimmed
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=fastqc_trimmed.%j

module load fastqc/0.11.9

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/

fastqc -t 12 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/*.fastq.gz

sbatch fastqc_polyg.sh 
Submitted batch job 9469286

while read file;
 do unzip $file 
done < filenames.txt

grep 'warn' */fastqc_data.txt > warnings.txt
grep 'fail' */fastqc_data.txt > failings.txt






grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_adaptertrimmed_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_adaptertrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' failings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_polygtrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' warnings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_polygtrimmed_fastas.txt

grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_finaltrim_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_finaltrim_fastas.txt


## Align and sort (make sorted bam files)
# first, get all finch names 
ls *_2P* | awk 'BEGIN {FS = "_"} {print $1}' | sort | uniq -c

for u in "${finch[@]}"; do
echo ${u}
bowtie2 -p 24 --phred33 --very-sensitive-local -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -I 149 -X 900 --rg-id ${u} --rg SM:${u} -1 ${u}.trimmed2_P.fq.gz -2 ${u}.trimmed2_2P.fq.gz -U ${u}.trimmed2_1U.fq.gz, ${u}.trimmed2_2U.fq.gz -S /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam 2>> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/alignment.output.txt
done


#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

module load bowtie2
module load picard
module load samtools
module list
 
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

echo "Job started on `date`"
# Download genome 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Camarhynchus_parvulus/latest_assembly_versions/GCF_901933205.1_STF_HiC/GCF_901933205.1_STF_HiC_genomic.fna.gz
#gzip -d GCF_901933205.1_STF_HiC_genomic.fna.gz 
#mv GCF_901933205.1_STF_HiC_genomic.fna GCF_901933205.fa
#bowtie2-build -f GCF_901933205.fa GCF_901933205.fa
#picard CreateSequenceDictionary R= GCF_901933205.fa O= GCF_901933205.fa.dict

#samtools faidx GCF_901933205.fa

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507 JP4481_round1 JP9655_round1 RHC507_round1 SM032_round1 SM059_round1 SM1067_round1 SM1200_round1 SM1240_round1 RHC097_round1 SM031_round1 SM040_round1 SM079_round1 SM1157_round1 SM1204_round1 SM1266_round1 JP5410_round1)

echo "Aligning fastas"


for u in "${finch[@]}"; do
echo ${u}
bowtie2 -p 24 --phred33 --very-sensitive-local -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -I 149 -X 900 --rg-id ${u} --rg SM:${u} -1 ${u}.trimmed2_1P.fq.gz -2 ${u}.trimmed2_2P.fq.gz -U ${u}.trimmed2_1U.fq.gz, ${u}.trimmed2_2U.fq.gz -S /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam 2>> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/alignment.output.txt
done

# convert to bams and sort
echo "Converting sams to bams"

for u in "${finch[@]}"; do 
samtools view -S -b /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam

echo "Sorting bams"

samtools sort /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam 
done


# some alignment stats 
echo "Computing stats on sorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

echo "Computing stats on unsorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/alignmentstats.txt 

echo "Job ended on `date`"


sbatch align_finaltrim.round1.sh 
Submitted batch job 9481231




#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.align.%j

module load bowtie2
module load picard
module load samtools
module list
 
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

echo "Job started on `date`"
# Download genome 
#wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Camarhynchus_parvulus/latest_assembly_versions/GCF_901933205.1_STF_HiC/GCF_901933205.1_STF_HiC_genomic.fna.gz
#gzip -d GCF_901933205.1_STF_HiC_genomic.fna.gz 
#mv GCF_901933205.1_STF_HiC_genomic.fna GCF_901933205.fa
#bowtie2-build -f GCF_901933205.fa GCF_901933205.fa
#picard CreateSequenceDictionary R= GCF_901933205.fa O= GCF_901933205.fa.dict

#samtools faidx GCF_901933205.fa

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)

echo "Aligning fastas"


for u in "${finch[@]}"; do
echo ${u}
bowtie2 -p 24 --phred33 --very-sensitive-local -x /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -I 149 -X 900 --rg-id ${u} --rg SM:${u} -1 ${u}_trimmed2_1P.fq.gz -2 ${u}_trimmed2_2P.fq.gz -U ${u}_trimmed2_1U.fq.gz, ${u}_trimmed2_2U.fq.gz -S /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam 2>> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/alignment.output.txt
done

# convert to bams and sort
echo "Converting sams to bams"

for u in "${finch[@]}"; do 
samtools view -S -b /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/samfiles/${u}.sam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam

echo "Sorting bams"

samtools sort /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam 
done


# some alignment stats 
echo "Computing stats on sorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

echo "Computing stats on unsorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/${u}.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/bamfiles/alignmentstats.txt 

echo "Job ended on `date`"

sbatch align_finaltrim.sh 
Submitted batch job 9483129


#!/bin/bash

#SBATCH --job-name=align_sequencess
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
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

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

finch=(JP9655 RHC507 SM031 SM032 SM040 SM059 RHC097 SM079 SM1204 SM1231 SM1237 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1157 SM1200 JP4481 SM1240 SM1266 SM1067 JP5410 SRR2917289 SRR2917290 SRR2917291 SRR2917292 SRR2917293 SRR2917294 SRR2917295 SRR2917296 SRR2917297 SRR2917298 SRR2917329 SRR2917330 SRR2917331 SRR2917332 SRR2917333 SRR2917334 SRR2917335 SRR2917336 SRR2917337 SRR2917338 SRR1607534 SRR1607532 SRR1607537 SRR1607539 SRR1607541 SRR1607533 SRR1607538 SRR1607540 SRR1607542 SRR1607535 SRR1607536 SRR1607504 SRR1607506 SRR1607505 SRR1607507)


echo "Computing stats on sorted bams"

for u in "${finch[@]}"; do
echo ${u}
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/${u}.sorted.bam
done > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedbamfiles/alignmentstats.txt 

sbatch sortedbam_stats.sh 
Submitted batch job 9515763



ls /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/

ls *all*bam | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt


## Clip overlapping read pairs
#!/bin/bash

#SBATCH --job-name=clipoverlap
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.clipping.%j

module load parallel

clipping () {
echo "clipping" "$@" >> /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/clippingstats.txt 

/home/u15/dannyjackson/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam --out /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.sorted.marked.clipped.bam --stats --params


echo "done " $@ >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/clippingstats.txt 
}


export -f clipping 

parallel -j 12 clipping :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt

sbatch clipping.all.sh 
Submitted batch job 1946694


## Compute statistics on bam files 
# some alignment stats 
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".depthstats

samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch stating_dupmarkedbams.sh 
Submitted batch job 1934432
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch stating_clipoverlap.sh 
Submitted batch job 9591145



cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/
ls *all* | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt

#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.sorted.marked.clipped.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.depthstats

samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt  





sbatch stating_all_dupmarkedbams.sh 
Submitted batch job 9544175


#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".all.sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".all.depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt  

sbatch stating_all_dupmarkedbams_fixerror.sh 
Submitted batch job 1934522


# compute average and sd depth per individual # sorted marked bam files
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 

sbatch averagedepth_sortedbams.sh 
Submitted batch job 9566373



# compute average and sd depth per individual # sorted marked bam files
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 



while read -r finch;
do 
  samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$finch".sorted.marked.clipped.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

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

while read -r finch;
do 
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 


# for all files 
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
  samtools index /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam

  apptainer exec ~/programs/gatk3_3.7-0.sif java -jar /usr/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa \
    -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
    -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

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
  -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".sorted.marked.clipped.bam \
  --targetIntervals /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelmaps/"$finch".intervals \
  -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

sbatch indelrealignment.sh 
Submitted batch job 1942426

# for all files 
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

sbatch ~/programs/slurmscripts/indelrealignment.all.sh

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

sbatch ~/programs/slurmscripts/indelmap.all.sh 
Submitted batch job 9724713

sbatch indelrealignment.all.sh 
Submitted batch job 9750831

# Compute per-sample sequencing depth for bases with mapping quality higher than 20


module load samtools
while read -r finch;
do 
  samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

ls *bam | awk 'BEGIN {FS = "."} {print $1}' | wc -l > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_indelrealigned.txt



### test for effects of duplicate sequencing effort










# python script for summarizing depth stats from ANGSD
python
import numpy as np
import pandas as pd 

df = pd.read_csv('stats.NC_044571.depthSample', sep='\t', header=None)
cols = df.shape[1]
multiplier = list(range(1, cols+1))

df_adj = df.mul(multiplier, axis = 1)

counts = df.sum(axis = 1)
totaldepth = df_adj.sum(axis = 1)
avgdepth = totaldepth / counts

df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = pd.concat([df_names, avgdepth], axis=1)

df_final = df_final.rename({0: 'avgdepth'}, axis=1)

df_final.groupby(['species', 'treatment']).mean()

df_final.sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'CRA'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'FOR'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'PAR'].sort_values(by=['avgdepth']).round(decimals=2)

print(df_final.groupby(['species', 'treatment']).size())

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')





CRA     post        6.712184
        pre        12.027188
FOR     post        5.676069
        pre         7.770516
PAR     post        6.038819
        pre        10.395772
