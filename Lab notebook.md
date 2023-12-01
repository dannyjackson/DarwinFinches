# Lab notebook

# November 6th
Goals: 1. Understand HPC storage and slurm processes, 2. Identify sequence names to species/treatment, 3. Write variant call script

Goal 1. Created HPC Notes.md and a local folder of HPC documentation. I feel really solid about my ability to use this and am ready to give it a try. I do want to talk to Sabrina about shared usage -- who else is running HPC things at the moment in the lab? If it's just me, I'm not going to worry too much about conserving time.

Goal 2. 

# these ones are all duplicate individuals:
SRR1607532	lamich_cra_PL15
SRR1607533	lamich_cra_PL15_b
SRR1607534	lamich_cra_PL16
SRR1607535	lamich_cra_PL16_b
SRR1607536	lamich_cra_PL16_c
SRR1607537	lamich_cra_PL4
SRR1607538	lamich_cra_PL4_b
SRR1607539	lamich_cra_PL7
SRR1607540	lamich_cra_PL7_b
SRR1607541	lamich_cra_PL9
SRR1607542	lamich_cra_PL9_b
SRR1607504	lamich_par_PARV1
SRR1607505	lamich_par_PARV1_b
SRR1607506	lamich_par_PARV2
SRR1607507	lamich_par_PARV2_b

# sabrina has already sorted and marked bam files
# there are 66, but only 59 in the google sheet
# it looks like that's because there are 7 additional files, which represent the sorted.marked bam files of the individuals with duplicate fastas
# i'd like to look through her code and see how she did this, but first get through the sanity check of it all
# create a subdirectory with all sorted.marked bam files but not the duplicates of single-duplicated individuals

cd /xdisk/mcnew/finch_wgs/fastqs
mkdir sorted_marked_bams
mv *sorted.marked.bam sorted_marked_bams
cd sorted_marked_bams

duplicatedindividuals=("SRR1607504.sorted.marked.bam" "SRR1607505.sorted.marked.bam" "SRR1607506.sorted.marked.bam" "SRR1607507.sorted.marked.bam" "SRR1607532.sorted.marked.bam" "SRR1607533.sorted.marked.bam" "SRR1607534.sorted.marked.bam" "SRR1607535.sorted.marked.bam" "SRR1607536.sorted.marked.bam" "SRR1607537.sorted.marked.bam" "SRR1607538.sorted.marked.bam" "SRR1607539.sorted.marked.bam" "SRR1607540.sorted.marked.bam" "SRR1607541.sorted.marked.bam" "SRR1607542.sorted.marked.bam")

for file in ${duplicatedindividuals[@]}; do
  mv $file ..
done

module load picard
module load samtools 
module load gatk

module load bcftools
module list


#!/bin/bash

#SBATCH --job-name=variant_calling
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.variant_calling.%j

module load bcftools

cd /xdisk/mcnew/finch_wgs/fastqs/danny

ref="/xdisk/mcnew/finch_wgs/fastqs/GCF_901933205.fa"
bamdir="/xdisk/mcnew/finch_wgs/fastqs/sorted_marked_bams/"
ID="darwinfinches"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.sorted.marked.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf


sbatch callvariants.slurm 
squeue --job 8529759

17284902

gatk HaplotypeCaller -R GCF_901933205.fa -I ${u}.sorted.marked.bam --native-pair-hmm-threads 10 -ERC GVCF -O ${u}.g.vcf >&log_${u} 

The other single individual genomes have this many lines:
75,512,578
185,176,268

cat /xdisk/mcnew/finch_wgs/fastqs/danny/darwinfinches_snps_multiallelic.vcf | wc -l
4,409,022

17 hours / 200
51,764,710 snps before it quits (not quick enough)
52,573,940


#!/bin/bash

#SBATCH --job-name=mergevcfs
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.mergevcfs.%j

module load bcftools

cd /xdisk/mcnew/finch_wgs/fastqs/danny

bcftools merge /xdisk/mcnew/finch_wgs/fastqs/*vcf.gz > combined.vcf

sbatch mergevcfs.slurm 
squeue --job 1679775


## November 8th
cat combined.vcf | wc -l
1,007,191,197

grep 'CHROM' combined.vcf
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	JP4481	JP5410	JP9655	lamich_PARV1	lamich_PARV2	lamich_PL15	lamich_PL16	lamich_PL4	lamich_PL7	lamich_PL9	RHC097	RHC507	SM031	SM032	SM040	SM059	SM079	SM1067	SM1083	SM1156	SM1157	SM1200	SM1204	SM1231	SM1237	SM1240	SM1266	SM1270	SM1271	SM1272	SM1273	SRR2917289	SRR2917290	SRR2917291	SRR2917292	SRR2917293	SRR2917294	SRR2917295	SRR2917296	SRR2917297	SRR2917298	SRR2917329	SRR2917330	SRR2917331	SRR2917332	SRR2917333	SRR2917334	SRR2917335	SRR2917336	SRR2917337	SRR2917338

cat darwinfinches_snps_multiallelic.vcf | wc -l
10,524,531

That's a wild number of snps. Is it because basically every site across the genome had at least one individual with missing data? Or because of something that's gone wrong in the bam generation etc?




#SBATCH --job-name=vcfstats
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.vcfstats.%j

module load plink

cd /xdisk/mcnew/finch_wgs/fastqs/danny/

plink --vcf /xdisk/mcnew/finch_wgs/fastqs/danny/combined.vcf --allow-extra-chr --missing --cluster-missing --freq

squeue --job 8537612
# the above one broke due to memory issues
squeue --job 8537759





# I'm thinking I'll maybe also redo the bam file process

cp /xdisk/mcnew/finch_wgs/fastqs/GCF_901933205.fa /xdisk/mcnew/dannyjackson/referencedata/

# copy this into a file of sample names

JP4481
JP5410
JP9655
lamich_PARV1
lamich_PARV2
lamich_PL15
lamich_PL16
lamich_PL4
lamich_PL7
lamich_PL9
RHC097
RHC507
SM031
SM032
SM040
SM059
SM079
SM1067
SM1083
SM1156
SM1157
SM1200
SM1204
SM1231
SM1237
SM1240
SM1266
SM1270
SM1271
SM1272
SM1273
SRR2917289
SRR2917290
SRR2917291
SRR2917292
SRR2917293
SRR2917294
SRR2917295
SRR2917296
SRR2917297
SRR2917298
SRR2917329
SRR2917330
SRR2917331
SRR2917332
SRR2917333
SRR2917334
SRR2917335
SRR2917336
SRR2917337
SRR2917338


# making the slurm script for running the align and sort script

#!/bin/bash

#SBATCH --job-name=alignsort
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.alignsort.%j

module load bwa
module load samtools
module load picard-tools

/home/u15/dannyjackson/programs/whole_genome_bioinformatics/align-and-sort.sh -t 12 -i /xdisk/mcnew/dannyjackson/reference_lists/samplenames.txt -p /xdisk/mcnew/dannyjackson/fastas/ -r /xdisk/mcnew/dannyjackson/reference_data/GCF_901933205.fa


SM1237_1.fastq.gz  SRR1607504_2.fastq.gz  SRR1607536_1.fastq.gz  SRR2917289_2.fastq.gz  SRR2917297_1.fastq.gz  SRR2917334_2.fastq.gz
JP9655_2.fastq.gz  SM1067_1.fastq.gz  SM1237_2.fastq.gz

sbatch ~/programs/slurmscripts/alignsort.slurm

squeue --job 8541861
I'm restarting it with 12 threads. Hoping this is much more time effective.

## November 9th



    alignsort <- function(sample) {

    samtools sort -T temp -@ 6 \
    -o "$sortedbamoutdir$SAMPLE"_sorted.bam \
    $bamoutdir$SAMPLE.bam

    picard-tools AddOrReplaceReadGroups \
    I=$sortedbamoutdir"$SAMPLE"_sorted.bam \
    O=$sortedbamoutdir"$SAMPLE"_sorted_RGadded.bam \
    RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=$SAMPLE

    picard-tools MarkDuplicates \
    I=$sortedbamoutdir"$SAMPLE"_sorted_RGadded.bam \
    O=$sortedbamoutdir"$SAMPLE"_sorted_RGadded_dupmarked.bam \
    M=$sortedbamoutdir"$SAMPLE".duplicate.metrics.txt

    samtools index \
    $sortedbamoutdir"$SAMPLE"_sorted_RGadded_dupmarked.bam

    rm $sortedbamoutdir"$SAMPLE"_sorted.bam
    rm $sortedbamoutdir"$SAMPLE"_sorted_RGadded.bam
    }


    mclapply("$seqs", align, mc.set.seed=FALSE)


sbatch ~/programs/slurmscripts/alignsort.slurm
8542479


module load bwa
module load samtools
module load picard

/home/u15/dannyjackson/programs/whole_genome_bioinformatics/align-and-sort.sh


#!/bin/bash

# shell script to assemble a bam file, then sort, mark duplicates, and index

./align.sh -t 12 -i /xdisk/mcnew/dannyjackson/reference_lists/samplenames.txt -p /xdisk/mcnew/dannyjackson/fastas/ -r /xdisk/mcnew/dannyjackson/reference_data/GCF_901933205.fa



module load parallel 
#!/bin/bash

# shell script to assemble a bam file, then sort, mark duplicates, and index

if [ $# -lt 1 ]
  then
    echo "Aligns fastq reads to a reference genome using bwa mem.
    Bam file is then sorted, duplicates are marked, and file is indexed using
    samtools and picard-tools.

    [-i] Sample list
    [-r] Reference genome

    OPTIONAL ARGUMENTS

    [-t] Number of threads to use
    [-p] Path to trimmed fastqs - the default is a directory called 'fastqs' as
         produced from the initial sorting
    [-b] Output directory for bam files - default is to make a directory
         called 'bam_files'
    [-s] Output directory for sorted bam files - default is to make a
         directory called 'sorted_bam_files'"

  else
    while getopts i:r:t:p:b:s: option
    do
    case "${option}"
    in
    i) seqs=${OPTARG};;
    r) ref=${OPTARG};;
    t) threads=${OPTARG};;
    p) fastqs_path=${OPTARG};;
    b) bamoutdir=${OPTARG};;
    s) sortbamoutdir=${OPTARG};;
    esac
    done

    threads="${threads:-1}"
    fastqs_path="${fastqs_path:-fastqs/}"
    bamoutdir="${bamoutdir:-bam_files/}"
    sortedbamoutdir="${sortedbamoutdir:-sorted_bam_files/}"
    if [ $bamoutdir == bam_files/ ]
      then
        mkdir bam_files
    fi

    echo "Beginning alignment for " $seqs>>bwa_alignment_log.txt

    align () {
    echo "aligning " $ID >> bwa_alignment_log.txt

    bwa mem -t $threads $ref $fastqs_path"$ID"_1.fastq.gz \
    $fastqs_path"$ID"_2.fastq.gz | \
    samtools view -b -o $bamoutdir$ID.bam -S 

    echo "sam file piped into samtools view to convert to .bam">>bwa_alignment_log.txt
    }
    
    export -f align 

    parallel align :::: "$seqs"
  fi


/home/u15/dannyjackson/programs/whole_genome_bioinformatics/align-and-sort.sh -t 12 -i /xdisk/mcnew/dannyjackson/reference_lists/samplenames.txt -p /xdisk/mcnew/dannyjackson/fastas/ -r /xdisk/mcnew/dannyjackson/reference_data/GCF_901933205.fa


## November 10th

variant call finished
20,447,414 lines

sbatch ~/programs/slurmscripts/vcfstats.slurm -i /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf

squeue --job 8546463


sbatch ~/programs/slurmscripts/filtervcf.slurm -i /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf.gz -n darwinfinches
squeue --job 8546467

sbatch ~/programs/slurmscripts/filtervcf.slurm -i /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf -n darwinfinches_nomaxDP
squeue --job 8599662

sed -i 's/\SRR2917/SRR/g' darwinfinches_filtered.recode.vcf

plink --vcf darwinfinches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out darwinfinches_filtered_mind2

# 9785604 variants and 51 people pass filters and QC.


plink --vcf darwinfinches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out darwinfinches_filtered_cleaned

cp darwinfinches_filtered_cleaned.vcf darwinfinches_filtered_cleaned_zip.vcf
bgzip darwinfinches_filtered_cleaned_zip.vcf

module load bcftools 

bcftools index darwinfinches_filtered_cleaned_zip.vcf.gz


# RAxML

~/programs/vcf2phylip/vcf2phylip.py -i darwinfinches_filtered_cleaned.vcf
mkdir pruned 

cd pruned 

module load python

python ~/programs/sula/filter_invariants_all.py ../darwinfinches_filtered_cleaned.min4.phy
mv variantsites.phy ../variantsites_mind2.phy 
mv variantsites_kept.txt ../variantsites_mind2_kept.txt 
cd .. 
rm -r pruned

cd /xdisk/mcnew/dannyjackson/finches/raxml

echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-1918469' > partitionfile.txt


#!/bin/bash

#SBATCH --job-name=raxml
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.raxml.%j

~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n darwinfinches -q /xdisk/mcnew/dannyjackson/finches/raxml/partitionfile.txt -s /xdisk/mcnew/dannyjackson/finches/vcfs/variantsites_mind2.phy -T 6 -p 12345 -N 10 ­-b 12345 -V

sbatch ~/programs/slurmscripts/raxml.slurm

squeue --job 8555921

# PCA

# library("devtools")
# install_github("zhengxwen/gdsfmt")
# install_github("zhengxwen/SNPRelate")

echo -e "CRA" > pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "CRA" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "FOR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 
echo -e "PAR" >> pops.txt 

#!/bin/bash

#SBATCH --job-name=PCA_all
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.PCA_all.%j

module load R 

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf  -o /xdisk/mcnew/dannyjackson/finches/PCA/all/ -p /xdisk/mcnew/dannyjackson/finches/PCA/all/pops.txt -n all -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/all/
sbatch ~/programs/slurmscripts/PCA_all.slurm
squeue --job 1682655

# subset by species
# cra
/xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf

#!/bin/bash

#SBATCH --job-name=PCA_cra
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.PCA_cra.%j

module load R 

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -o /xdisk/mcnew/dannyjackson/finches/PCA/cra/ -p /xdisk/mcnew/dannyjackson/finches/PCA/cra/pops.txt -n cra -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/cra/
sbatch ~/programs/slurmscripts/PCA_cra.slurm
squeue --job 1682677

post
post
pre
pre
pre
pre
pre
post
post
post
post
post
post

JP4481  CRA post
JP5410  CRA post
lamich_PL15 CRA pre
lamich_PL16 CRA pre
lamich_PL4  CRA pre
lamich_PL7  CRA pre
lamich_PL9  CRA pre
SM1067  CRA post
SM1156  CRA post
SM1157  CRA post
SM1200  CRA post
SM1240  CRA post
SM1266  CRA post

# par
bcftools view -s JP9655,lamich_PARV1,lamich_PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf

#!/bin/bash

#SBATCH --job-name=PCA_cra
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.PCA_cra.%j

module load R 

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf  -o /xdisk/mcnew/dannyjackson/finches/PCA/par/ -p /xdisk/mcnew/dannyjackson/finches/PCA/par/pops.txt -n par -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/par/
sbatch ~/programs/slurmscripts/PCA_par.slurm
squeue --job 1682674

post
pre
pre
post
post
post
post
post
post
post
pre
pre
pre
pre
pre
pre
pre
pre
pre
pre


JP9655  PAR post
lamich_PARV1  PAR pre
lamich_PARV2  PAR pre
RHC097  PAR post
RHC507  PAR post
SM031 PAR post
SM032 PAR post
SM040 PAR post
SM059 PAR post
SM079 PAR post
SRR329  PAR pre
SRR330  PAR pre
SRR331  PAR pre
SRR332  PAR pre
SRR333  PAR pre
SRR334  PAR pre
SRR335  PAR pre
SRR336  PAR pre
SRR337  PAR pre
SRR338  PAR pre

# for
bcftools view -s SM1083,SM1204,SM1231,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf



post
post
post
post
post
post
post
post
pre
pre
pre
pre
pre
pre
pre
pre
pre
pre

#!/bin/bash

#SBATCH --job-name=PCA_for
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.PCA_for.%j

module load R 

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf  -o /xdisk/mcnew/dannyjackson/finches/PCA/for/ -p /xdisk/mcnew/dannyjackson/finches/PCA/for/pops.txt -n for -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/for/
sbatch ~/programs/slurmscripts/PCA_for.slurm

squeue --job 1682737

SM1083  FOR post
SM1204  FOR post
SM1231  FOR post
SM1237  FOR post
SM1270  FOR post
SM1271  FOR post
SM1272  FOR post
SM1273  FOR post
SRR289  FOR pre
SRR290  FOR pre
SRR291  FOR pre
SRR292  FOR pre
SRR293  FOR pre
SRR294  FOR pre
SRR295  FOR pre
SRR296  FOR pre
SRR297  FOR pre
SRR298  FOR pre

















# subset by pre/post

# pre
bcftools view -s lamich_PARV1,lamich_PARV2,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/pre.vcf

# post

bcftools view -s JP4481,JP5410,JP9655,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SM1067,SM1083,SM1156,SM1157,SM1200,SM1204,SM1231,SM1237,SM1240,SM1266,SM1270,SM1271,SM1272,SM1273 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/post.vcf


# pre sample list 
PAR
PAR
CRA
CRA
CRA
CRA
CRA
FOR
FOR
FOR
FOR
FOR
FOR
FOR
FOR
FOR
FOR
PAR
PAR
PAR
PAR
PAR
PAR
PAR
PAR
PAR
PAR

# post sample list
CRA
CRA
PAR
PAR
PAR
PAR
PAR
PAR
PAR
PAR
CRA
FOR
CRA
CRA
CRA
FOR
FOR
FOR
CRA
CRA
FOR
FOR
FOR
FOR

# full sample list 
PAR = Camarhynchus parvulus	Small tree finch
CRA = Platyspiza crassirostris	Vegetarian finch
FOR = Geospiza fortis	Medium ground finch

JP4481  CRA post
JP5410  CRA post
JP9655  PAR post
lamich_PARV1  PAR pre
lamich_PARV2  PAR pre
lamich_PL15 CRA pre
lamich_PL16 CRA pre
lamich_PL4  CRA pre
lamich_PL7  CRA pre
lamich_PL9  CRA pre
RHC097  PAR post
RHC507  PAR post
SM031 PAR post
SM032 PAR post
SM040 PAR post
SM059 PAR post
SM079 PAR post
SM1067  CRA post
SM1083  FOR post
SM1156  CRA post
SM1157  CRA post
SM1200  CRA post
SM1204  FOR post
SM1231  FOR post
SM1237  FOR post
SM1240  CRA post
SM1266  CRA post
SM1270  FOR post
SM1271  FOR post
SM1272  FOR post
SM1273  FOR post
SRR289  FOR pre
SRR290  FOR pre
SRR291  FOR pre
SRR292  FOR pre
SRR293  FOR pre
SRR294  FOR pre
SRR295  FOR pre
SRR296  FOR pre
SRR297  FOR pre
SRR298  FOR pre
SRR329  PAR pre
SRR330  PAR pre
SRR331  PAR pre
SRR332  PAR pre
SRR333  PAR pre
SRR334  PAR pre
SRR335  PAR pre
SRR336  PAR pre
SRR337  PAR pre
SRR338  PAR pre



# PCA output 
JP = Jeff Potost (umass amherst)



# these are all labelled FOR but one is clustering with CRA
SM1231


# this one is labelled CRA but is clustering with FOR
SM1156



# Wrongly labelled samples (maybe)

SM1156 (labelled CRA, clusters with FOR) > DEFINITELY FORTIS
  got swapped with 1231

SM1231 (labelled FOR, is actually CRA)

either SM1270 or SRR293 is odd in the PAR

# in each species specific PCA, there's one odd individual:
# CRA: one wierd; it's SM1156, is actually FOR
# FOR: one weird; it's SM1231, is actually CRA
# PAR: one weird; it's JP9655, is actually CRA






## November 14th

I wasn't sure I wanted to reassign individuals before I got the raxml results, but I really want to see an fst plot...


# revised full sample list 
PAR = Camarhynchus parvulus	Small tree finch
CRA = Platyspiza crassirostris	Vegetarian finch
FOR = Geospiza fortis	Medium ground finch

JP4481  CRA post
JP5410  CRA post
JP9655  CRA post
lamich_PARV1  PAR pre
lamich_PARV2  PAR pre
lamich_PL15 CRA pre
lamich_PL16 CRA pre
lamich_PL4  CRA pre
lamich_PL7  CRA pre
lamich_PL9  CRA pre
RHC097  PAR post
RHC507  PAR post
SM031 PAR post
SM032 PAR post
SM040 PAR post
SM059 PAR post
SM079 PAR post
SM1067  CRA post
SM1083  FOR post
SM1156  FOR post
SM1157  CRA post
SM1200  CRA post
SM1204  FOR post
SM1231  CRA post
SM1237  FOR post
SM1240  CRA post
SM1266  CRA post
SM1270  FOR post
SM1271  FOR post
SM1272  FOR post
SM1273  FOR post
SRR289  FOR pre
SRR290  FOR pre
SRR291  FOR pre
SRR292  FOR pre
SRR293  FOR pre
SRR294  FOR pre
SRR295  FOR pre
SRR296  FOR pre
SRR297  FOR pre
SRR298  FOR pre
SRR329  PAR pre
SRR330  PAR pre
SRR331  PAR pre
SRR332  PAR pre
SRR333  PAR pre
SRR334  PAR pre
SRR335  PAR pre
SRR336  PAR pre
SRR337  PAR pre
SRR338  PAR pre


# PCA All
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "CRA" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "FOR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt
echo -e "PAR" >> pops.txt



#!/bin/bash

#SBATCH --job-name=PCA_all
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.PCA_all.%j

module load R 

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf  -o /xdisk/mcnew/dannyjackson/finches/PCA/all/ -p /xdisk/mcnew/dannyjackson/finches/PCA/all/pops.txt -n all -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/all/
sbatch ~/programs/slurmscripts/PCA_all.slurm
squeue --job 8561923




# FST

grep 'NC_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf | awk '{print $1}' | sort -u | wc -l

NC_044571.1
NC_044572.1
NC_044573.1
NC_044574.1
NC_044575.1
NC_044576.1
NC_044577.1
NC_044578.1
NC_044579.1
NC_044580.1
NC_044581.1
NC_044582.1
NC_044583.1
NC_044584.1
NC_044585.1
NC_044586.1
NC_044587.1
NC_044588.1
NC_044589.1
NC_044590.1
NC_044591.1
NC_044592.1
NC_044593.1
NC_044594.1
NC_044595.1
NC_044596.1
NC_044597.1
NC_044598.1
NC_044599.1
NC_044600.1
NC_044601.1

grep 'NW_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf | awk '{print $1}' | sort -u | wc -l

grep -e 'NW_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.chroms.vcf

# PAR

cd /xdisk/mcnew/dannyjackson/finches/fst/par

echo -e "lamich_PARV1" > pre.txt 
echo -e "lamich_PARV2" >> pre.txt 
echo -e "RHC097" >> post.txt 
echo -e "RHC507" >> post.txt 
echo -e "SM031" >> post.txt 
echo -e "SM032" >> post.txt 
echo -e "SM040" >> post.txt 
echo -e "SM059" >> post.txt 
echo -e "SM079" >> post.txt 
echo -e "SRR329" >> pre.txt 
echo -e "SRR330" >> pre.txt 
echo -e "SRR331" >> pre.txt 
echo -e "SRR332" >> pre.txt 
echo -e "SRR333" >> pre.txt 
echo -e "SRR334" >> pre.txt 
echo -e "SRR335" >> pre.txt 
echo -e "SRR336" >> pre.txt 
echo -e "SRR337" >> pre.txt 
echo -e "SRR338" >> pre.txt 



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

sed -i 's/NC_//g' out.windowed.weir.fst
sed -i 's/NW_//g' out.windowed.weir.fst

library(qqman)
fst<-read.table("out.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "par_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt 

head -1 out.weir.fst > chroms.weir.fst
grep 'NC' out.weir.fst >> chroms.weir.fst

sed -i 's/NC_//g' chroms.weir.fst

sed -i 's/NC_//g' out.weir.fst
sed -i 's/NW_//g' out.weir.fst


library(qqman)
fst<-read.table("chroms.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "par.chroms.fst.snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()




# CRA
cd /xdisk/mcnew/dannyjackson/finches/fst/cra

echo -e "JP4481" > post.txt 
echo -e "JP5410" >> post.txt 
echo -e "JP9655" >> post.txt 
echo -e "lamich_PL15" > pre.txt 
echo -e "lamich_PL16" >> pre.txt 
echo -e "lamich_PL4" >> pre.txt 
echo -e "lamich_PL7" >> pre.txt 
echo -e "lamich_PL9" >> pre.txt 
echo -e "SM1067" >> post.txt 
echo -e "SM1157" >> post.txt 
echo -e "SM1200" >> post.txt 
echo -e "SM1231" >> post.txt 
echo -e "SM1240" >> post.txt 
echo -e "SM1266" >> post.txt 


module load vcftools
module load R


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

head -1 out.windowed.weir.fst > chroms.windowed.weir.fst
grep 'NC' out.windowed.weir.fst >> chroms.windowed.weir.fst

sed -i 's/NC_//g' chroms.windowed.weir.fst

library(qqman)
fst<-read.table("chroms.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "cra_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt 

head -1 out.weir.fst > chroms.weir.fst
grep 'NC' out.weir.fst >> chroms.weir.fst

sed -i 's/NC_//g' chroms.weir.fst


library(qqman)
fst<-read.table("chroms.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "cra.chroms.fst.snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()


# FOR
cd /xdisk/mcnew/dannyjackson/finches/fst/for

echo -e "SM1083" > post.txt 
echo -e "SM1156" >> post.txt 
echo -e "SM1204" >> post.txt 
echo -e "SM1237" >> post.txt 
echo -e "SM1270" >> post.txt 
echo -e "SM1271" >> post.txt 
echo -e "SM1272" >> post.txt 
echo -e "SM1273" >> post.txt 
echo -e "SRR289" > pre.txt 
echo -e "SRR290" >> pre.txt 
echo -e "SRR291" >> pre.txt 
echo -e "SRR292" >> pre.txt 
echo -e "SRR293" >> pre.txt 
echo -e "SRR294" >> pre.txt 
echo -e "SRR295" >> pre.txt 
echo -e "SRR296" >> pre.txt 
echo -e "SRR297" >> pre.txt 
echo -e "SRR298" >> pre.txt 




module load vcftools
module load R


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

head -1 out.windowed.weir.fst > chroms.windowed.weir.fst
grep 'NC' out.windowed.weir.fst >> chroms.windowed.weir.fst

sed -i 's/NC_//g' chroms.windowed.weir.fst

library(qqman)
fst<-read.table("chroms.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]

xu <- mean(fstsubset$WEIGHTED_FST)
s <- sd(fstsubset$WEIGHTED_FST)
fstsubset$ZFST = (fstsubset$WEIGHTED_FST - xu)/s

SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "for_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()

pdf(file = "for_1000window_zfst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="ZFST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt 

head -1 out.weir.fst > chroms.weir.fst
grep 'NC' out.weir.fst >> chroms.weir.fst

sed -i 's/NC_//g' chroms.weir.fst


library(qqman)
fst<-read.table("chroms.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "for.chroms.fst.snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()







## Thursday November 16th

# maybe refilter each VCF just by species?

# par
bcftools view -s lamich_PARV1,lamich_PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf

bcftools view -s RHC097,RHC507,SM031,SM032,SM040,SM059,SM079 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.vcf

bcftools view -s lamich_PARV1,lamich_PARV2,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.vcf



# cra
bcftools view -s JP4481,JP5410,JP9655,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf

bcftools view -s lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf


bcftools view -s JP4481,JP5410,JP9655,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf


# for
bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf

bcftools view -s SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.vcf

bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.vcf


echo -e "SM1083" > post.txt 
echo -e "SM1156" >> post.txt 
echo -e "SM1204" >> post.txt 
echo -e "SM1237" >> post.txt 
echo -e "SM1270" >> post.txt 
echo -e "SM1271" >> post.txt 
echo -e "SM1272" >> post.txt 
echo -e "SM1273" >> post.txt 
echo -e "SRR289" > pre.txt 
echo -e "SRR290" >> pre.txt 
echo -e "SRR291" >> pre.txt 
echo -e "SRR292" >> pre.txt 
echo -e "SRR293" >> pre.txt 
echo -e "SRR294" >> pre.txt 
echo -e "SRR295" >> pre.txt 
echo -e "SRR296" >> pre.txt 
echo -e "SRR297" >> pre.txt 
echo -e "SRR298" >> pre.txt 



Sort out transformed FST plots


n = dim(Q(obj.snmf, K = 2))[1]
fst.values[fst.values<0] = 0.000001
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))



# I think at this point it's useful to compare pre- post- statistics
Coverage:
Nucleotide diversity:


plink --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf  --allow-extra-chr --missing --cluster-missing --freq


plink --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf  --cluster-missing
plink --file mydata --test-missing



# Friday November 17th or Monday November whatevertheheck
https://www.york.ac.uk/res/dasmahapatra/teaching/MBiol_sequence_analysis/workshop4_2019.html#nucleotide_diversity



RAxML finished running -- species cluster without bias by pre-post group. Geospiza fortis did not form a monophyletic clade and instead contained parvulus. I didn't use an outgroup so this analysis is somewhat limited.




# Tuesday November 28th
PCA finished running
Double check that CRA are properly assigned pre-post. Two aren't clustering. Looks like either SM1200 or JP9655 and PL9. I think they're wrong -- PL9 looks like it's assigned "post" but i think it's gotta be a pre.


vcftools --vcf L_donovani_all_samples.vcf --window-pi 10000 --out L_donovani_all_samples_10kb




# Same week, various notes up through Thursday November 30th
I really want one clean and clear script that I can run that will run genomewide fst, nucleotide diversity, tajima's d, sweed, and bayescan and create one final file where each row is a 



# cra 
# pre
cd /xdisk/mcnew/dannyjackson/finches/nucleotidediversity/cra/pre

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf --window-pi 10000

head -1 out.windowed.pi > chroms.windowed.pi
grep 'NC' out.windowed.pi >> chroms.windowed.pi

sed -i 's/NC_//g' chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' chroms.windowed.pi | column -t > chroms.windowed.pi.formanhattan




R

pi.all <- read.table("chroms.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = "cra_pre_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "cra_pre_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

# this gives regions of significance 

mydf[ which(mydf$myBg=='TRUE'),]

write.csv(mydf[ which(mydf$myBg=='TRUE'),], "cra_pre_td_sig.csv", row.names=FALSE)

SNP myBg   CHROM BIN_START BIN_END N_VARIANTS        PI
94646 94646 TRUE 44601.1   4400001 4410000        307 0.0161473

grep '44601.1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene'



# post
cd /xdisk/mcnew/dannyjackson/finches/nucleotidediversity/cra/post

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf --window-pi 10000

head -1 out.windowed.pi > chroms.windowed.pi
grep 'NC' out.windowed.pi >> chroms.windowed.pi

sed -i 's/NC_//g' chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' chroms.windowed.pi | column -t > chroms.windowed.pi.formanhattan




R

pi.all <- read.table("chroms.windowed.pi",header=T)
pi.subset<-pi.all[complete.cases(pi.all),]

SNP<-c(1: (nrow(pi.subset)))

lower = min(pi.subset$PI)
upper = max(pi.subset$PI)
cutoff = upper - ((upper-lower)*0.05)

LessThanCutoff <- pi.subset$PI < cutoff

myBg <- !LessThanCutoff



mydf<-data.frame(SNP,myBg,pi.subset)

pdf(file = "cra_post_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(pi.subset$PI,br=20)
dev.off()

pdf(file = "cra_post_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(PI ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()

# this gives regions of significance 

mydf[ which(mydf$myBg=='TRUE'),]


SNP myBg   CHROM BIN_START BIN_END N_VARIANTS        PI
70708 70708 TRUE 44583.1   2200001 2210000        371 0.0162050
95064 95064 TRUE 44601.1   4400001 4410000        311 0.0160899

grep '44583.1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene'

# this is the only gene that isn't 
NC_044583.1     Gnomon  gene    2203974 2217325 .       -       .       ID=gene-LOC115908533;Dbxref=GeneID:115908533;Name=LOC115908533;gbkey=Gene;gene=LOC115908533;gene_biotype=protein_coding

protocadherin gamma-A12-like

grep '44601.1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene'



./selection_scans.sh -n cra -o /xdisk/mcnew/dannyjackson/finches/nucleotidediversity/cra/ -p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf -q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf

    n) name=${OPTARG};;
    o) outDir=${OPTARG};;
    p) pop1=${OPTARG};;
    q) pop2=${OPTARG};;
 


# i need to search GFF for genes in the chromosome of column "CHROM" of the CSV (col3)

# testing

# this finds shared values between file1$1 and file2$2

awk -F"[,\t]" 'NR==FNR{a[$1]=$1","$2; next} ($2 in a){print a[$2]","$1}' file1.txt file2.txt


# i need to find shared values between file1$3 and file2$1

# remake test files
1,Brian,1000
4,Jason,1010
8,Nick,400
13,Sean,410

1000 3044	
400  4466
1010 1206

awk -F"[,\t]" 'NR==FNR{a[$3]=$3","$1; next} ($1 in a){print a[$1]","$3}' file1.txt file2.txt

# so that works, but it's printing a combination of file2$1 and file1$3
# what I need is the entire line of file 2 if it matches file1$3

awk -F"[,\t]" 'NR==FNR{a[$3]=$3","$1; next} ($1 in a){print $0}' test2.test.csv file2.txt

awk -F"[,\t]" 'NR==FNR{a[$3]=$3","$1; next} ($1 in a){print $0}' file1.txt file2.txt


# generate test csv 
awk NR\>1 cra.pi_sig.csv > test.csv
awk 'BEGIN {FS = ",";OFS = ","} $3="NC_0"$3' test.csv > test2.csv


# generate test giff

awk '$0 !~ /\#/' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | head > test.gff
grep '44583.1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | grep 'LOC115908533' >> test.gff

# modify awk script to apply to my test files
awk -F"[,\t]" 'NR==FNR{a[$3]=$3","$1; next} ($1 in a){print $0}' test2.csv test.gff

# okay this works now! try it with the non-test files

awk -F"[,\t]" 'NR==FNR{a[$3]=$3","$1; next} ($1 in a){print $0}' test2.csv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff

# that worked! now I have to get it to identify not just where file1$3 is a match but also where file1$4 < file2$5 and file1$5 > file2$4


awk '$5>= 2217320 && $5>= 2203970' test.gff

# the above works! I need to incorporate it into the other awk script now ugh


awk -F"[,\t]" 'NR==FNR{a[$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' test2.csv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' 

# that works too! I now, annoyingly, need to incorporate the modifications that I made to the output csv

awk -F"[,\t]" 'NR==FNR{a[NC_0$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' test2.csv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' 

## biiiiiiiiiiitch yes!!!!!!!

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' cra.pi_sig.csv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' 


awk NR\>1 cra.pi_sig.csv > test.csv
awk 'BEGIN {FS = ",";OFS = ","} $3="NC_0"$3' test.csv > test2.csv



# tomorrow, continue to build out the script incorporating pre and post pops into the nucleotide diversity scans, then do fst, tajima's d, sweed, and bayescan. The last one is gonna be the most annoying i think. The final step should be to incorporate all of these output into one single file that has rows of genes with columns that attribute which programs identified them.




awk -F"[,\t]" 'NR==FNR{a[NC_0$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}'  /xdisk/mcnew/dannyjackson/finches/nucleotidediversity/cra/interestpop/cra.pi_sig.csv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' > ${outDir}/interestpop/${name}.pi_sig_genes.csv


# Friday December 1st
/xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf




## previous fst scripts
# PAR

module load vcftools
module load R

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

head -1 out.windowed.weir.fst > chroms.windowed.weir.fst
grep 'NC' out.windowed.weir.fst >> chroms.windowed.weir.fst

sed -i 's/NC_//g' chroms.windowed.weir.fst



library(qqman)
fst<-read.table("chroms.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]

xu <- mean(fstsubset$WEIGHTED_FST)
s <- sd(fstsubset$WEIGHTED_FST)
fstsubset$ZFST = (fstsubset$WEIGHTED_FST - xu)/s

SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "par_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()

pdf(file = "par_1000window_zfst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="ZFST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()

# CRA

module load vcftools
module load R


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

head -1 out.windowed.weir.fst > chroms.windowed.weir.fst
grep 'NC' out.windowed.weir.fst >> chroms.windowed.weir.fst

sed -i 's/NC_//g' chroms.windowed.weir.fst


library(qqman)
fst<-read.table("chroms.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]

xu <- mean(fstsubset$WEIGHTED_FST)
s <- sd(fstsubset$WEIGHTED_FST)
fstsubset$ZFST = (fstsubset$WEIGHTED_FST - xu)/s

SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "cra_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()

pdf(file = "cra_1000window_zfst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="ZFST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()


# FOR
module load vcftools
module load R


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt --fst-window-size 1000 --fst-window-step 1000 

head -1 out.windowed.weir.fst > chroms.windowed.weir.fst
grep 'NC' out.windowed.weir.fst >> chroms.windowed.weir.fst

sed -i 's/NC_//g' chroms.windowed.weir.fst

library(qqman)
fst<-read.table("chroms.windowed.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]

xu <- mean(fstsubset$WEIGHTED_FST)
s <- sd(fstsubset$WEIGHTED_FST)
fstsubset$ZFST = (fstsubset$WEIGHTED_FST - xu)/s

SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "for_1000window.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()

pdf(file = "for_1000window_zfst.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="BIN_START",p="ZFST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst", cex = 0.2))
dev.off()



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf --weir-fst-pop pre.txt --weir-fst-pop post.txt 

head -1 out.weir.fst > chroms.weir.fst
grep 'NC' out.weir.fst >> chroms.weir.fst

sed -i 's/NC_//g' chroms.weir.fst


library(qqman)
fst<-read.table("chroms.weir.fst", header=TRUE)
fstsubset<-fst[complete.cases(fst),]
SNP<-c(1: (nrow(fstsubset)))
mydf<-data.frame(SNP,fstsubset)

pdf(file = "for.chroms.fst.snps.pdf", width = 20, height = 7, useDingbats=FALSE)
print(manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="SNP",logp=FALSE,ylab="Weighted Weir and Cockerham Fst"))
dev.off()




#!/bin/bash

#SBATCH --job-name=selectionscans
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.selectionscans.%j

module load vcftools
module load R

chmod +x ~/programs/DarwinFinches/selection_scans.sh

~/programs/DarwinFinches/selection_scans.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -n cra -o /xdisk/mcnew/dannyjackson/finches/cra/ -p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf -q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf -r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_pops.txt -s /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_pops.txt -g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff

chmod -x ~/programs/DarwinFinches/selection_scans.sh




# Tajima's D

vcftools --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_pyrr.filtered.geno25.maf1.vcf --TajimaD 20000 


sed -i 's/VYXE//g' out.Tajima.D

awk '{sub(/\./,"",$1)}1' out.Tajima.D | column -t > out.Tajima.D.formanhattan


R

taj.all <- read.table("out.Tajima.D",header=T)
taj.subset<-taj.all[complete.cases(taj.all),]

SNP<-c(1: (nrow(taj.subset)))

lower = min(taj.subset$TajimaD)
upper = max(taj.subset$TajimaD)
lower_cutoff = lower + ((upper-lower)*0.025)
upper_cutoff = upper - ((upper-lower)*0.025)

MoreThanLower <- taj.subset$TajimaD > lower_cutoff
LessThanUpper <- taj.subset$TajimaD < upper_cutoff
significant <- MoreThanLower & LessThanUpper 

myBg <- !significant


mydf<-data.frame(SNP,myBg,taj.subset)

sigdf <-  mydf[which(mydf$myBg),]

write.table(sigdf, file = "sigtd.tsv")

pdf(file = "pyrr_urban_td_hist.pdf", width = 10, height = 5, useDingbats=FALSE)
hist(taj.subset$TajimaD,br=20)
dev.off()

pdf(file = "pyrr_urban_td.pdf", width = 20, height = 7, useDingbats=FALSE)

plot(TajimaD ~ SNP, col= "white", pch = 21, bg=ifelse(myBg  == "TRUE", 'red', 'gray'),
     data = mydf,
     xaxt = "n", bty = "l", xlab = "chr", cex = 1)
# add custom axis labels
axis(1, at = mydf$SNP,
     labels = mydf$CHROM,
     las = 2)

dev.off()




# Next steps:
Consider https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-07095-8#Sec8
which used zfst, zHp, pi-ratio, and also:
XP-EHH: Cross-population extended haplotype homozygosity

Rubin et al. defined and applied a Z-score test for heterozygosity depression (ZHp) on chicken genome sequence, which basically expresses how much the expected heterozygosity in chromosome windows deviate from the average genome heterozygosity [11].

https://www.nature.com/articles/nature08832

Run ARGwaver (Ancestral Recombination Graph)


# ideal header for output of gene table
genename,chromosome,startpos,endpos,cra,par,for,fst,pi,etc

python --version
# If you have trouble installing pixy in an environment using python 3.9, try rolling back to python 3.8.



conda config --add channels conda-forge
conda install -c conda-forge pixy
conda install -c bioconda htslib
# test installation
pixy --help
