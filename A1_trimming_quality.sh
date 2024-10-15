# trimming and quality assessment

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

