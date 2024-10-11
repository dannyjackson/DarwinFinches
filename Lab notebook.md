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

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par.autosomes.vcf

bcftools view -s RHC097,RHC507,SM031,SM032,SM040,SM059,SM079 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.autosomes.vcf

bcftools view -s lamich_PARV1,lamich_PARV2,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.autosomes.vcf


# cra

bcftools view -s JP4481,JP5410,JP9655,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra.autosomes.vcf

bcftools view -s lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.autosomes.vcf

bcftools view -s JP4481,JP5410,JP9655,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.autosomes.vcf


# for
bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for.autosomes.vcf


bcftools view -s SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.autosomes.vcf

bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.vcf

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.autosomes.vcf

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

~/programs/DarwinFinches/selection_scans.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -n cra -o /xdisk/mcnew/dannyjackson/finches/cra/ -p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf -q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf -r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_pops.txt -s /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_pops.txt -g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff - /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_bayescan_popfile.txt 

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


conda create --name pixy
conda activate pixy
conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib
pip install numpy==1.22.1
# test installation
pixy --help


Consider paring back the script so that it doesn't have to do everything / have a dozen options passed to it. Idk just a thought.

## December 6th
Goals:
1. Confirm that I used the trimmed sequences in my pipeline. If not, start it from the beginning with the proper ones.
2. Then work on the script to run various sliding window summary stats 

Retrimming, starting all analyses from scratch. Good thing I am making my scripts repeatable instead of doing all this shit by hand <3 

awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' cra.tD_sig.tsv /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff  | grep 'ID=gene' 


awk -F"[,\t]" '{print "NC_0"substr($3, 1, length($2)-1)".1" }' 

awk -F"[ \t]" '{print $4}' cra.tD_sig.tsv 

substr($3, 1, length($2)-1)


> ${outDir}/tajimasd/interestpop/${name}.tD_sig_genes.csv

echo 'IND,POP' > /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_popgenome_popfile.txt
cat /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_bayescan_popfile.txt >> /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_popgenome_popfile.txt
# sorting out dxy
# https://markravinet.github.io/Chapter8.html
library(tidyverse)
# library(devtools)
# install_github("pievos101/PopGenome")
library(PopGenome)

cra <- readData("/xdisk/mcnew/dannyjackson/finches/cra/dxy", format = "VCF", include.unknown = TRUE, FAST = TRUE)

get.sum.data(cra)
cra@n.biallelic.sites 
# 3050107
cra@populations

cra_info <- read_delim("/xdisk/mcnew/dannyjackson/finches/reference_lists/cra_popgenome_popfile.txt", delim = ",")

populations <- split(cra_info$IND, cra_info$POP)

cra <- set.populations(cra, populations, diploid = T)

cra@populations


cra_sw <- sliding.window.transform(cra, width = 100000, jump = 25000)

cra_sw <- diversity.stats(cra_sw, pi = TRUE)
cra_sw <- F_ST.stats(cra_sw, mode = "nucleotide")

nd <- cra_sw@nuc.diversity.within/100000
pops <- c("post", "pre")
# set population names
colnames(nd) <- paste0(pops, "_pi")

fst <- t(cra_sw@nuc.F_ST.pairwise)
dxy <- get.diversity(cra_sw, between = T)[[2]]/100000

# get column names 
x <- colnames(fst)
# replace all occurrences of pop1 with house
# x <- sub("pop1", "post", x)
# does the same thing as above but by indexing the pops vector
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)

colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

cra_data <- as.tibble(data.frame(nd, fst, dxy))


~/programs/DarwinFinches/selection_scans.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -n cra -o /xdisk/mcnew/dannyjackson/finches/cra/ -p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf -q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf -r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_pops.txt -s /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_pops.txt -g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff - /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_bayescan_popfile.txt 



# Got pixy installed!

sed 's/\,/\t/g' /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_bayescan_popfile.txt > /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pixy_popfile.txt

module load samtools
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf.gz



conda activate pixy

module load samtools

pixy --stats fst,dxy --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf.gz --populations /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pixy_popfile.txt --window_size 1000 --output_folder /xdisk/mcnew/dannyjackson/finches/cra/dxy --bypass_invariant_check yes

sed -i 's/NC_//g' ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi

sbatch ~/programs/slurmscripts/pixy_cra.slurm 
Submitted batch job 1706575


awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/interestpop/${name}.tD_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/interestpop/${name}.tD_sig_genes.csv


/xdisk/mcnew/dannyjackson/finches/cra/dxy



## Friday December 8th
1. Revise selection scan script to use pixy for pi, dxy, and fst

Stuff from previous script that I'm removing:

# reference population analysis
vcftools --vcf ${pop1} --window-pi 10000 --out ${outDir}/nucleotidediversity/referencepop/${name}

head -1 ${outDir}/nucleotidediversity/referencepop/${name}.windowed.pi > ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi
grep 'NC' ${outDir}/nucleotidediversity/referencepop/${name}.windowed.pi >> ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi | column -t > ${outDir}/nucleotidediversity/referencepop/${name}.chroms.windowed.pi.formanhattan

Rscript ~/programs/DarwinFinches/nucleotidediversity.r ${outDir}/nucleotidediversity/referencepop ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/nucleotidediversity/referencepop/${name}.pi_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/nucleotidediversity/referencepop/${name}.pi_sig_genes.csv

# population of interest analysis
vcftools --vcf ${pop2} --window-pi 10000 --out ${outDir}/nucleotidediversity/interestpop/${name}

head -1 ${outDir}/nucleotidediversity/interestpop/${name}.windowed.pi > ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi
grep 'NC' ${outDir}/nucleotidediversity/interestpop/${name}.windowed.pi >> ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi

sed -i 's/NC_//g' ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi

awk '{sub(/\./,"",$1)}1' ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi | column -t > ${outDir}/nucleotidediversity/interestpop/${name}.chroms.windowed.pi.formanhattan

Rscript ~/programs/DarwinFinches/nucleotidediversity.r ${outDir}/nucleotidediversity/interestpop ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$3]=$0; b=$4; c=$5; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/nucleotidediversity/interestpop/${name}.pi_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/nucleotidediversity/interestpop/${name}.pi_sig_genes.csv



# fst 
vcftools --vcf ${vcf} --weir-fst-pop ${p1file} --weir-fst-pop ${p2file} --out ${outDir}/fst/${name} --fst-window-size 10000 


head -1 ${outDir}/fst/${name}.windowed.weir.fst > ${outDir}/fst/${name}.chroms.windowed.weir.fst
grep 'NC' ${outDir}/fst/${name}.windowed.weir.fst >> ${outDir}/fst/${name}.chroms.windowed.weir.fst

sed -i 's/NC_//g' ${outDir}/fst/${name}.chroms.windowed.weir.fst

Rscript ~/programs/DarwinFinches/fstscans.r ${outDir}/fst ${name}

awk -F"[,\t]" 'NR==FNR{a["NC_0"$2]=$0; b=$3; c=$4; next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/fst/${name}.zfst_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/fst/${name}.zfst_sig_genes.csv






## working notes for the revision
~/programs/DarwinFinches/selection_scans.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -n cra -o /xdisk/mcnew/dannyjackson/finches/cra/ -p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf -q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf -r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_pops.txt -s /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_pops.txt -g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff 



awk -F"[ \t]" 'NR==FNR{a["NC_0"substr($4, 1, length($4)-1)".1"]=$0; b=$5; c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/tajimasd/interestpop/${name}.tD_sig.csv ${gff} | grep 'ID=gene' > ${outDir}/tajimasd/interestpop/${name}.tD_sig_genes.csv


head -1 pixy_dxy.txt > cra.chroms.pixy_dxy.txt
grep 'NC' pixy_dxy.txt >> cra.chroms.pixy_dxy.txt


awk '{sub(/\./,"",$1)}1' cra.chroms.pixy_dxy.txt | column -t > cra.chroms.pixy_dxy.formanhattan.txt

sbatch ~/programs/DarwinFinches/pixy.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf.gz -o /xdisk/mcnew/dannyjackson/finches/cra/ -r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pixy_popfile.txt -w 10000


sbatch ~/programs/DarwinFinches/pixy.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf.gz -o /xdisk/mcnew/dannyjackson/finches/par/ -r /xdisk/mcnew/dannyjackson/finches/reference_lists/par_pixy_popfile.txt -w 10000

sbatch ~/programs/DarwinFinches/pixy.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf.gz -o /xdisk/mcnew/dannyjackson/finches/for/ -r /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pixy_popfile.txt -w 10000

outDir="/xdisk/mcnew/dannyjackson/finches/cra"
name="cra"

head -1 ${outDir}/pixy/pixy_dxy.txt > ${outDir}/pixy/${name}.chroms.pixy_dxy.txt
grep 'NC' ${outDir}/pixy/pixy_dxy.txt >> ${outDir}/pixy/${name}.chroms.pixy_dxy.txt
sed -i 's/NC_//g' ${outDir}/pixy/${name}.chroms.pixy_dxy.txt
awk '{sub(/\./,"",$1)}1' ${outDir}/pixy/${name}.chroms.pixy_dxy.txt | column -t > ${outDir}/pixy/${name}.chroms.pixy_dxy.formanhattan.txt


Rscript ~/programs/DarwinFinches/pixy.r ${outDir}/pixy ${name}


## December 12th

awk -F"[,\t]" 'NR==FNR{a["NC_0"$4]=$0; b=($5-10000); c=($5 +20000); next} ($1 in a && $5 >= b && $4<=c){print $0}' ${outDir}/pixy/fst/${name}.zfst_sig.csv ${gff} | grep 'ID=gene'


# CRA
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.vcf.gz

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.vcf.gz \
-n cra \
-o /xdisk/mcnew/dannyjackson/finches/cra/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1709717 ## look in /vcfs for output

# FOR
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.vcf.gz

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.vcf.gz \
-n for \
-o /xdisk/mcnew/dannyjackson/finches/for/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.vcf \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.vcf \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1709719

# PAR
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.vcf.gz

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.vcf.gz \
-n par \
-o /xdisk/mcnew/dannyjackson/finches/par/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.vcf \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.vcf \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/par_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1711788








# fixing PCA script

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf -o /xdisk/mcnew/dannyjackson/finches/PCA/cra/ -p /xdisk/mcnew/dannyjackson/finches/PCA/cra/pops.txt -n cra -s y


dataset <- "/xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf"
outDir <- "/xdisk/mcnew/dannyjackson/finches/PCA/cra/"
pops <- "/xdisk/mcnew/dannyjackson/finches/PCA/cra/pops.txt"
name <- "cra"

library(gdsfmt)
library(SNPRelate)

#open the GDS file
genofile <- snpgdsOpen(paste0(outDir,"/",name,".gds"))

#get population information
pop_code <- scan(paste(pops), what=character())
table(pop_code)

#prune SNPS based on linkage disequilibrium
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, autosome.only = FALSE)
snpset.id <- unlist(snpset)

#run PCA
#with subsetted snps
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = FALSE)

#to calculate percent of variation that is accounted for by top PCA components:
pc.percent <- pca$varprop*100
sink(paste0(outDir,"/",name,"_percentvariation.txt"), append=TRUE, split=FALSE)
head(round(pc.percent, 2))
sink()

#get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

#make data frame
tab <- data.frame(sample.id = pca$sample.id,
pop = factor(pop_code)[match(pca$sample.id, sample.id)],
EV1 = pca$eigenvect[,1],    # the first eigenvector
EV2 = pca$eigenvect[,2],    # the second eigenvector
stringsAsFactors = FALSE)

#save files:
write.table(tab, file = paste0(outDir,"/",name,"_subset.txt"), sep = "\t")

tab <- read.table(paste0(outDir,"/",name,"_subset.txt"), header = TRUE, sep = "\t", stringsAsFactors=T)

PC1.PV = pc.percent[1]
PC2.PV = pc.percent[2]

#draw it:
pdf(file = paste0(outDir,"/",name,"_pca_populations_subsetted.pdf"), useDingbats=FALSE)

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), xlab=paste0("PC1 (Percent Variation =",PC1.PV,")"), ylab=paste0("PC2 (Percent Variation =",PC1.PV,")"))
legend("topright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))

dev.off()

pdf(file = paste0(outDir,"/",name,"_pca_individuals_subsetted.pdf"), useDingbats=FALSE)

plot(tab$EV1, tab$EV2, col=as.integer(tab$sample.id), xlab="PC1", ylab="PC2")
legend("topright", legend=levels(tab$sample.id), pch="o", col=1:nlevels(tab$sample.id))

dev.off()





## December 13th 

I need to make sense of the strange patterns on the Z chromosome. The first step is to identify which individuals are which sex. The next is to redo all previous analyses using only autosomal chromosomes (after presenting at lab meeting today)

module load bcftools
module load samtools

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.vcf

tabix /xdisk/mcnew/dannyjackson/fin/darwinfinches_filtered.geno25.maf1.notindep.vcf.gz

bcftools view -r NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.zchrom.vcf


module load plink

plink --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.zchrom.vcf --allow-extra-chr --missing --cluster-missing --freq

bcftools stats -d /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.zchrom.vcf

module load vcftools 

vcftools --het --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.zchrom.vcf

vcftools --het --gzvcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.vcf.gz

# F should be around 0 for males and 1 for females on the Z chromosome

awk 'FNR==NR{a[FNR]=$5;next};{$NF=a[FNR]};1' genomic.het zchrom.het


awk 'NR==FNR {a[FNR]=$5; next} {print a[FNR] "\t" $5; }' genomic.het zchrom.het > bothF.het

sed -i 1d bothF.het

awk '{print $2 / $1}' bothF.het


bcftools stats /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.vcf.gz


plink --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.vcf.gz --allow-extra-chr --missing --freq








~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.notindep.zchrom.vcf -o /xdisk/mcnew/dannyjackson/finches/zchrom/PCA -p /xdisk/mcnew/dannyjackson/finches/PCA/all/pops.txt -n all_zchrom -s y



## To Do After Lab Meeting
Follow up on filtering and depth
Why isn't my vcf retaining depth scores?
Calculate depth by individual and depth by Z score





# January 19th


# CRA
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.notindep.autosomes.vcf



bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.autosomes.vcf.gz 

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/cra.notindep.autosomes.vcf.gz \
-n cra \
-o /xdisk/mcnew/dannyjackson/finches/cra/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.notindep.autosomes.vcf.gz \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.notindep.autosomes.vcf.gz \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1773800 ## look in ~ for output

# FOR
bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.autosomes.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.autosomes.vcf.gz

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.autosomes.vcf.gz

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/for.notindep.autosomes.vcf.gz \
-n for \
-o /xdisk/mcnew/dannyjackson/finches/for/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.notindep.autosomes.vcf.gz \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.notindep.autosomes.vcf.gz \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1775262 ## look in ~ for output

# PAR
bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.autosomes.vcf.gz

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.vcf.gz

bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.vcf.gz > /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.autosomes.vcf

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.autosomes.vcf.gz

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.autosomes.vcf
tabix /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.autosomes.vcf.gz 

sbatch ~/programs/DarwinFinches/selection_scans.sh \
-v /xdisk/mcnew/dannyjackson/finches/vcfs/par.notindep.autosomes.vcf.gz \
-n par \
-o /xdisk/mcnew/dannyjackson/finches/par/ \
-p /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.notindep.autosomes.vcf.gz \
-q /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.notindep.autosomes.vcf.gz \
-r /xdisk/mcnew/dannyjackson/finches/reference_lists/par_pixy_popfile.txt \
-g /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff \
-w 10000

Submitted batch job 1775483 ## look in ~ for output


# February 15th
bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf  > darwinfinches_qualitysort.vcf

bcftools view  -i  'MIN(FORMAT/DP)>5'  /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.vcf  > darwinfinches_DP5.vcf

bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf > b10k_filtered.recode.stats.txt

bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > b10k_filtered.geno25.maf1.stats.txt

vcftools --vcf /xdisk/mcnew/dannyjackson/vcfs/ --min-meanDP 2 --remove-indels --recode --out b10k_filtered

sed -i 's/\SRR2917/SRR/g' darwinfinches_filtered.recode.vcf

module load picard
module load samtools 
module load gatk

module load bcftools
module load vcftools

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --min-meanDP 4 --remove-indels --recode --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth4.vcf
# After filtering, kept 13030803 out of a possible 13808013 Sites

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --min-meanDP 6 --remove-indels --recode --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth6.vcf
# After filtering, kept 9715863 out of a possible 13808013 Sites

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --min-meanDP 8 --remove-indels --recode --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth8.vcf

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --min-meanDP 10 --remove-indels --recode --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth10.vcf
# After filtering, kept 313275 out of a possible 13808013 Sites

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --min-meanDP 12 --remove-indels --recode --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth12.vcf
# After filtering, kept 190501 out of a possible 13808013 Sites


# February 16th 

bcftools stats -v /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth6.vcf.recode.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth6.stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth6.vcf.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysortdepth6.depthstats

vcftools --vcf darwinfinches_snps_multiallelic.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic


darwinfinches_snps_multiallelic.idepth
R
data <- read.csv('darwinfinches_snps_multiallelic.idepth', sep ='\t')
samplelist <- read.csv('samplelist.txt', sep ='\t')
> mean(data$MEAN_DEPTH)
[1] 6.454152

df <- merge(data,samplelist,by="INDV") 

library("dplyr")                             
 

 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))
 
print(group_mean)

  TREATMENT Mean_Depth
  <chr>          <dbl>
1 post            4.13
2 pre            8.52


group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(SPECIES,TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(Mean_Depth = mean))

print(group_mean)


  SPECIES TREATMENT Mean_Depth
  <chr>   <chr>          <dbl>
1 CRA     post            4.05
2 CRA     pre            10.4 
3 FOR     post            4.66
4 FOR     pre             6.85
5 PAR     post            3.63
6 PAR     pre             9.15



# now looking at ones where we filtered on read depth:
darwinfinches_qualitysortdepth6.depthstats.idepth

library("dplyr")                             

data <- read.csv('darwinfinches_qualitysortdepth6.depthstats.idepth', sep ='\t')
samplelist <- read.csv('samplelist.txt', sep ='\t')
df <- merge(data,samplelist,by="INDV") 

group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))

print(group_mean)

1 post            4.61
2 pre             9.82

group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(SPECIES,TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(Mean_Depth = mean))

print(group_mean)

  SPECIES TREATMENT Mean_Depth
  <chr>   <chr>          <dbl>
1 CRA     post            4.52
2 CRA     pre            12.1 
3 FOR     post            5.18
4 FOR     pre             7.89
5 PAR     post            4.07
6 PAR     pre            10.5 


# try filtering on depth with bcftools

bcftools view  -i  'MIN(FMT/DP)>6' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort_depth6_bcftools.vcf

bcftools view  -i  'MIN(FMT/DP)>6 & MIN(FMT/GQ)>100'  darwinfinches_snps_multiallelic.vcf > darwinfinches_qualitydepth


bcftools stats darwinfinches_snps_multiallelic.vcf > darwinfinches_snps_multiallelic.stats
bcftools stats darwinfinches_qualitysort.vcf > darwinfinches_qualitysort.stats

vcftools --vcf darwinfinches_qualitysort.vcf --depth 
13808013 sites

vcftools --vcf darwinfinches_snps_multiallelic.vcf --depth 
20446037 sites

bgzip /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf

bcftools index /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf.gz

bcftools view  -i  'MIN(FMT/DP)>6'  /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf.gz  > darwinfinches_quality_DP6.vcf



bcftools view  -i 'FMT/AD[GT] > 6' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf.gz  > darwinfinches_quality_DP6_alt.vcf

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.vcf --minDP 2 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP2.vcf
20446037

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP2.vcf.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP2.depth

vcftools --vcf darwinfinches_vcftools.DP2.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP2.miss


vcftools --vcf darwinfinches_snps_multiallelic.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP4.vcf

vcftools --vcf darwinfinches_vcftools.DP4.vcf.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP4.depth 

vcftools --vcf darwinfinches_vcftools.DP4.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP4.miss



R
data.raw <- read.csv('darwinfinches_snps_multiallelic.imiss', sep ='\t')
data.dp2 <- read.csv('darwinfinches_vcftools.DP2.miss.imiss', sep ='\t')
data.dp4 <- read.csv('darwinfinches_vcftools.DP4.miss.imiss', sep ='\t')



samplelist <- read.csv('samplelist.txt', sep ='\t')

mean(data.raw$F_MISS)
0.03721052

mean(data.dp2$F_MISS)
0.09974047

mean(data.dp4$F_MISS)
0.3071487


library("dplyr")   

df.raw <- merge(data.raw,samplelist,by="INDV") 
                          
 group_mean <- df.raw %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))

group_mean

1 post      0.0497
2 pre       0.0261

df.dp2 <- merge(data.dp2,samplelist,by="INDV") 
                          
 group_mean <- df.dp2 %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))

group_mean

1 post      0.150 
2 pre       0.0547

df.dp4 <- merge(data.dp4,samplelist,by="INDV") 
                          
 group_mean <- df.dp4 %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))

group_mean

1 post       0.482
2 pre        0.152

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.vcf --minDP 10 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_vcftools.DP10.vcf


vcftools --vcf darwinfinches_DP2.vcf --depth
176660

vcftools --vcf darwinfinches_AD2.vcf --depth
25800

vcftools --vcf darwinfinches_DP5.vcf --depth
103838

vcftools --vcf darwinfinches_DP5_AD5.vcf --depth
13938 sites

vcftools --vcf darwinfinches_AD5.vcf --depth
13938 sites

interactive -a mcnew -n 12 -t 3-00:00 #bayescan

ref="/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
bamdir="/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/sorted_marked_bams/"
ID="darwinfinches"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.sorted.marked.bam | bcftools call -mv -V indels > "$ID"_snps_multiallelic.vcf

# revised full sample list 
PAR = Camarhynchus parvulus	Small tree finch
CRA = Platyspiza crassirostris	Vegetarian finch
FOR = Geospiza fortis	Medium ground finch

INDV	SPECIES	TREATMENT
JP4481	CRA	post
JP5410	CRA	post
JP9655	CRA	post
lamich_PARV1	PAR	pre
lamich_PARV2	PAR	pre
lamich_PL15	CRA	pre
lamich_PL16	CRA	pre
lamich_PL4	CRA	pre
lamich_PL7	CRA	pre
lamich_PL9	CRA	pre
RHC097	PAR	post
RHC507	PAR	post
SM031	PAR	post
SM032	PAR	post
SM040	PAR	post
SM059	PAR	post
SM079	PAR	post
SM1067	CRA	post
SM1083	FOR	post
SM1156	FOR	post
SM1157	CRA	post
SM1200	CRA	post
SM1204	FOR	post
SM1231	CRA	post
SM1237	FOR	post
SM1240	CRA	post
SM1266	CRA	post
SM1270	FOR	post
SM1271	FOR	post
SM1272	FOR	post
SM1273	FOR	post
SRR2917289	FOR	pre
SRR2917290	FOR	pre
SRR2917291	FOR	pre
SRR2917292	FOR	pre
SRR2917293	FOR	pre
SRR2917294	FOR	pre
SRR2917295	FOR	pre
SRR2917296	FOR	pre
SRR2917297	FOR	pre
SRR2917298	FOR	pre
SRR2917329	PAR	pre
SRR2917330	PAR	pre
SRR2917331	PAR	pre
SRR2917332	PAR	pre
SRR2917333	PAR	pre
SRR2917334	PAR	pre
SRR2917335	PAR	pre
SRR2917336	PAR	pre
SRR2917337	PAR	pre
SRR2917338	PAR	pre








# notes to work on next


# sweed 

bcftools view -s 'NOCA003,NOCA004,NOCA006,NOCA008,NOCA012,NOCA013' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf.gz > /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf

plink --vcf /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid  --out urban_noca.filtered.geno25.maf1

cd northerncardinals

~/programs/sweed/SweeD -name urban_noca -input /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf -grid 100 -length 100000


bgzip /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf
tabix /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz
zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_urban.txt

scp -r to_b10k/omega dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed

grep 'Position' SweeD_Report.noca_urban | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_urban_manhattan


while IFS=$'\t' read -r col1 col2 col3 col4 col5
do 
    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
        scaf="${col1:2}"
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$scaf $col1 $col2 $col3 $col4 $col5" >> SweeD_Report.noca_urban_manhattan
    fi
done < SweeD_Report.urban_noca

grep 'Position' SweeD_Report.urban_noca | head -n 1 > SweeD_Report.header
awk '{print $0, "Scaffold"}' SweeD_Report.header | awk '{print $6,$1,$2,$3,$4,$5}' > SweeD_Report.noca_urban_manhattan_scaffnames

while read -r col1 col2 col3 col4 col5 col6;
do 
 if [[ $col1 = 'Scaffold' ]]; then  
        continue
  else
    scaff=$(awk "FNR == ${col1} {print}" scaff_names_noca_urban.txt)
    echo "$scaff $col2 $col3 $col4 $col5 $col6" >> SweeD_Report.noca_urban_manhattan_scaffnames
  fi
done < SweeD_Report.noca_urban_manhattan

sed -i 's/VYXE//g' SweeD_Report.noca_urban_manhattan_scaffnames
# sed -i 's/\(.\{1\}\).1/\1/' SweeD_Report.noca_urban_manhattan_scaffnames
# awk '{sub(/\./,"",$1)}1' SweeD_Report.noca_urban_manhattan_scaffnames | column -t > SweeD_Report.noca_urban_manhattan_scaffnames.formanhattan

R
library(qqman)
sweed<-read.table("SweeD_Report.noca_urban_manhattan_scaffnames", header=TRUE)
sweed.subset<-sweed[complete.cases(sweed),]
SNP<-c(1: (nrow(sweed.subset)))

lower = min(sweed.subset$Likelihood)
upper = max(sweed.subset$Likelihood)
cutoff = upper - ((upper-lower)*0.05)
LessThanCutoff <- sweed.subset$Likelihood < cutoff
myBg <- !LessThanCutoff
mydf<-data.frame(SNP,myBg,sweed.subset)
sigdf <-  mydf[which(mydf$myBg),]
write.table(sigdf, file = "sigsweed.tsv")

mydf<-data.frame(SNP,sweedsubset)

pdf(file = "noca_urban_sweed.pdf", width = 20, height = 7, useDingbats=FALSE)
    print(manhattan(mydf,chr="Scaffold",bp="Position",p="Likelihood",snp="Position",logp=FALSE,ylab="CLR"))
dev.off()

scp dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/cardinalis/to_b10k/sweed/northerncardinals/noca_urban_sweed.pdf .





# February 19th 
vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_qualitysort.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.vcf



vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.vcf.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing6.vcf

# 12489107 sites out of a possible 13808013 Sites
# 0 is more permissive, 1 is less permissive. This means that we included sites up to 0.4 freq missingness, I believe. 

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.vcf.recode.vcf --max-missing 0.45 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing45.vcf

# 13134533 out of a possible 13808013 Sites


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.vcf.recode.vcf --max-missing 1 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing1.vcf

# 107610

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.vcf.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing8.vcf

# 2330074 out of a possible 13808013

vcftools --vcf darwinfinches.quality.DP4.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4


vcftools --vcf darwinfinches.quality.DP4.missing6.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing45.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing45

vcftools --vcf darwinfinches.quality.DP4.missing8.vcf.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing8

vcftools --vcf darwinfinches.quality.DP4.missing8.vcf.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches.quality.DP4.missing8

darwinfinches.quality.DP4.missing8.idepth


R
data.dp4.qual.miss45 <- read.csv('darwinfinches.quality.DP4.missing45.imiss', sep ='\t')

data.dp4.qual.miss6 <- read.csv('darwinfinches.quality.DP4.missing6.imiss', sep ='\t')
data.dp4.qual.miss8 <- read.csv('darwinfinches.quality.DP4.missing8.imiss', sep ='\t')
data.dp4.qual.miss8.depth <- read.csv('darwinfinches.quality.DP4.missing8.idepth', sep ='\t')

data.dp4.qual <- read.csv('darwinfinches.quality.DP4.imiss', sep ='\t')


samplelist <- read.csv('samplelist.txt', sep ='\t')

mean(data.dp4.qual$F_MISS)
0.2922971

mean(data.dp4.qual.miss8$F_MISS)
0.1625411

mean(data.dp4.qual.miss6$F_MISS)
0.25847

mean(data.dp4.qual.miss45$F_MISS)
0.2683867

mean(data.raw$F_MISS)
0.03721052

mean(data.dp2$F_MISS)
0.09974047

mean(data.dp4$F_MISS)
0.3071487


mean(data.dp4.qual.miss8.depth$MEAN_DEPTH)
# 9.827197

library("dplyr")   

df <- merge(data.dp4.qual.miss8.depth,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))

group_mean
1 post            6.57
2 pre            12.7 


 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(SPECIES,TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))

group_mean

1 CRA     post            6.48
2 CRA     pre            15.0 
3 FOR     post            7.32
4 FOR     pre            10.7 
5 PAR     post            5.82
6 PAR     pre            13.5 


df <- merge(data.dp4.qual.miss8,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))

group_mean

1 post      0.290 
2 pre       0.0494


 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(SPECIES,TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))
group_mean

1 CRA     post      0.294 
2 CRA     pre       0.0115
3 FOR     post      0.231 
4 FOR     pre       0.0991
5 PAR     post      0.351 
6 PAR     pre       0.0237

df <- merge(data.dp4.qual.miss6,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))

group_mean

1 post      0.444 
2 pre       0.0935


group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(SPECIES,TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(Mean_Depth = mean))

print(group_mean)



# February 20th
Goals for today:
Split raw vcf into species
Filter and calculate stats
Meet with Sabrina about the filtering step and make a decision about what to move forward with


sed 's/\SRR2917/SRR/g' darwinfinches_snps_multiallelic.vcf > darwinfinches_snps_multiallelic.renamed.vcf

# for
bcftools view -s SM1156,SM1083,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.renamed.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.vcf


bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.vcf

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing8


vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing6.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing6.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing6

# 11,967,809 sites
head -1 /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt > /xdisk/mcnew/dannyjackson/finches/vcfs/for/samplelist.txt
grep 'FOR' /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt >> /xdisk/mcnew/dannyjackson/finches/vcfs/for/samplelist.txt

sed -i 's/\SRR2917/SRR/g' /xdisk/mcnew/dannyjackson/finches/vcfs/for/samplelist.txt


R
data.miss <- read.csv('for.quality.DP4.missing6.imiss', sep ='\t')
data.depth <- read.csv('for.quality.DP4.missing6.idepth', sep ='\t')
samplelist <- read.csv('/xdisk/mcnew/dannyjackson/finches/vcfs/for/samplelist.txt', sep ='\t')

mean(data.miss$F_MISS)
# 0.2463031
mean(data.depth$MEAN_DEPTH)
# 6.407272

library("dplyr")   

df <- merge(data.miss,samplelist,by="INDV") 
            
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))
group_mean

# post       0.353
# pre        0.161

df <- merge(data.depth,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))
group_mean

# post            5.03
# pre             7.51


# par
bcftools view -s lamich_PARV1,lamich_PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.renamed.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.vcf


bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.vcf

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.missing6.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.missing6.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/par/par.quality.DP4.missing6

# 12,473,991 sites

head -1 /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt > /xdisk/mcnew/dannyjackson/finches/vcfs/par/samplelist.txt
grep 'PAR' /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt >> /xdisk/mcnew/dannyjackson/finches/vcfs/par/samplelist.txt

sed -i 's/\SRR2917/SRR/g' /xdisk/mcnew/dannyjackson/finches/vcfs/par/samplelist.txt


R
data.miss <- read.csv('par.quality.DP4.missing6.imiss', sep ='\t')
data.depth <- read.csv('par.quality.DP4.missing6.idepth', sep ='\t')
samplelist <- read.csv('/xdisk/mcnew/dannyjackson/finches/vcfs/par/samplelist.txt', sep ='\t')

mean(data.miss$F_MISS)
# 0.2240445
mean(data.depth$MEAN_DEPTH)
# 7.702896

library("dplyr")   

df <- merge(data.miss,samplelist,by="INDV") 
            
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))
group_mean

# post      0.518 
# pre       0.0526

df <- merge(data.depth,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))
group_mean

# post            3.90
# pre             9.92


# cra
bcftools view -s JP4481,JP5410,JP9655,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.renamed.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.vcf

bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.vcf

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.missing6.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.missing6

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.missing6.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/cra/cra.quality.DP4.missing6

# 9,997,683 out of a possible 13808013 Sites


head -1 /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt > /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt
grep 'CRA' /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt >> /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt

sed -i 's/\SRR2917/SRR/g' /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt


R
data.miss <- read.csv('cra.quality.DP4.missing6.imiss', sep ='\t')
data.depth <- read.csv('cra.quality.DP4.missing6.idepth', sep ='\t')
samplelist <- read.csv('/xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt', sep ='\t')

mean(data.miss$F_MISS)
# 0.2593543

mean(data.depth$MEAN_DEPTH)
# 7.123363

library("dplyr")   

df <- merge(data.miss,samplelist,by="INDV") 
            
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))
group_mean

# post      0.392 
# pre       0.0213

df <- merge(data.depth,samplelist,by="INDV") 
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))
group_mean

# post            4.64
# pre            11.6 

# revised full sample list 
PAR = Camarhynchus parvulus	Small tree finch
lamich_PARV1	PAR	pre
lamich_PARV2	PAR	pre
RHC097	PAR	post
RHC507	PAR	post
SM031	PAR	post
SM032	PAR	post
SM040	PAR	post
SM059	PAR	post
SM079	PAR	post
SRR2917329	PAR	pre
SRR2917330	PAR	pre
SRR2917331	PAR	pre
SRR2917332	PAR	pre
SRR2917333	PAR	pre
SRR2917334	PAR	pre
SRR2917335	PAR	pre
SRR2917336	PAR	pre
SRR2917337	PAR	pre
SRR2917338	PAR	pre

CRA = Platyspiza crassirostris	Vegetarian finch
JP4481	CRA	post
JP5410	CRA	post
JP9655	CRA	post
lamich_PL15	CRA	pre
lamich_PL16	CRA	pre
lamich_PL4	CRA	pre
lamich_PL7	CRA	pre
lamich_PL9	CRA	pre
SM1067	CRA	post
SM1157	CRA	post
SM1200	CRA	post
SM1231	CRA	post
SM1240	CRA	post
SM1266	CRA	post

FOR = Geospiza fortis	Medium ground finch
SM1270	FOR	post
SM1271	FOR	post
SM1272	FOR	post
SM1273	FOR	post
SRR2917289	FOR	pre
SRR2917290	FOR	pre
SRR2917291	FOR	pre
SRR2917292	FOR	pre
SRR2917293	FOR	pre
SRR2917294	FOR	pre
SRR2917295	FOR	pre
SRR2917296	FOR	pre
SRR2917297	FOR	pre
SRR2917298	FOR	pre
SM1083	FOR	post
SM1156	FOR	post
SM1204	FOR	post
SM1237	FOR	post

INDV	SPECIES	TREATMENT





FOR
# 18 individuals
Worst case scenarios:
0.6: (18*0.6 = 10.8) 9 pre and 2 post
0.75: (18*0.75 = ) 9 pre and 6 post
0.8: (18*0.8 = 14.4) 9 pre and 6 post

# 11,967,809 sites with 60%
# 3,291,439 sites with 80%

Missingness
# 0.2463031

Depth
# 6.407272

Missingness
# post       0.353
# pre        0.161

Depth
# post            5.03
# pre             7.51


PAR
# 12,473,991 sites
# 19 individuals, 12 pre and 7 post 
# (potentially drop the 3 )
Worst case scenarios:
0.6: (19*0.6 = 11.4) 12 pre and 0 post
0.8: (19*0.8 = 15.2) 12 pre and 4 post

Missingness
# 0.2240445

Depth
# 7.702896

Frequency
# post      0.518 
# pre       0.0526

Depth
# post            3.90
# pre             9.92


CRA
# 14 individuals, 5 pre and 9 post
Worst case scenarios:
0.6: (14*0.6 = 8.4) 5 pre and 4 post
0.6: (14*0.6 = 8.4) 9 post and 0 pre

0.8: (14*0.8 = 11.2) 5 pre and 8 post


# 9,997,683 out of a possible 13808013 Sites

Missingness
# 0.2593543

Depth
# 7.123363

Missingness
# post      0.392 
# pre       0.0213

Depth
# post            4.64
# pre            11.6 



# Next steps:
1. Incorporate first low coverage sequence data for post-samples from Sabrina
2. Split by species
3. Filter by quality >100 
     # justin uses a score of > 20
4. Filter by depth of 4
5. Use 80% missingness as our filtering step
Keep notes on how many snps get kept in each step.
Repeat analyses with depth of 3 and 60% missingness.



# February 23rd
1. Incorporate first low coverage sequence data for post-samples from Sabrina

sed -i 's/.fq.gz//g' filenames.txt

while read file;
do
mv ${file}.fq.gz ${file}.round1.fq.gz
done < filenames.txt


sbatch adapterremoval.sh 
Submitted batch job 9183050

sbatch adapterremoval.2.sh  
Submitted batch job 9183067


sbatch align.sh
Submitted batch job 9190364

sbatch markduplicates.sh 
Submitted batch job 9222374

sbatch markduplicates_sm059.sh 
Submitted batch job 9249610


sbatch variantcalling.slurm 
Submitted batch job 9250235

# DADI notes
module load anaconda/2022.05

# installation

pip install dadi
pip install -U numpy==1.22.0

# make pop file 
SM1156    post
SM1083    post
SM1204    post
SM1237    post
SM1270    post
SM1271    post
SM1272    post
SM1273    post
SRR289    pre
SRR290    pre
SRR291    pre
SRR292    pre
SRR293    pre
SRR294    pre
SRR295    pre
SRR296    pre
SRR297    pre
SRR298    pre


python
import dadi

dd = dadi.Misc.make_data_dict_vcf("/xdisk/mcnew/dannyjackson/finches/vcfs/for/for.quality.DP4.missing8.recode.vcf", "for_popfile.txt")
fs = dadi.Spectrum.from_data_dict(dd, ['pre', 'post'], projections = [20, 30], polarized = False)

scrambled = fs.scramble()




# March 5th

### ANGSD
module load htslib/1.19.1

git clone https://github.com/ANGSD/angsd.git 
cd angsd 
make

# quick start 
ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/sorted_marked_bams/ > /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

sed -i 's/^/\/xdisk\/mcnew\/dannyjackson\/finch_wgs\/danny_fastqs\/sorted_marked_bams\//g' /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

while read -r bam;
do
    samtools index "$bam" &
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt
wait

# calculate depth 
#!/bin/bash

#SBATCH --job-name=angsd_depth
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_depth.%j

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -doDepth 1 -out depth.alls 1 -nThreads 12

sbatch angsd_depth.sh 
Submitted batch job 9308935


# make beagle 

#!/bin/bash

#SBATCH --job-name=angsd_makebeagle
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=500gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_makebeagle.%j

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/darwinfinches -nThreads 10 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 1e-6 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt


sbatch angsd_makebeagle.sh 
Submitted batch job 9308948

# filtering options
-setMinDepth [int] -minMapQ [int] -minQ [int] 
~/programs/angsd/angsd  -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -doDepth 1 -out depth.all -doCounts 1 -minMapQ 100 -minQ 100













## depth after including all reads

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_snps_multiallelic_allreads.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_snps_multiallelic_allreads

# 9,997,683 out of a possible 13808013 Sites


head -1 /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt > /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt
grep 'CRA' /xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt >> /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt

sed -i 's/\SRR2917/SRR/g' /xdisk/mcnew/dannyjackson/finches/vcfs/cra/samplelist.txt


R
data.depth <- read.csv('darwinfinches_snps_multiallelic_allreads.idepth', sep ='\t')
samplelist <- read.csv('/xdisk/mcnew/dannyjackson/finches/vcfs/samplelist.txt', sep ='\t')


mean(data.depth$MEAN_DEPTH)
# 4.779342

library("dplyr")   

df <- merge(data.miss,samplelist,by="INDV") 
            
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(F_MISS),
                 list(F_MISS = mean))
group_mean

# post
# pre

df <- merge(data.depth,samplelist,by="INDV")
                          
 group_mean <- df %>%
    # Specify group indicator, column, function
    group_by(TREATMENT) %>%
    # Calculate the mean of the "Frequency" column for each group
    summarise_at(vars(MEAN_DEPTH),
                 list(MEAN_DEPTH = mean))
group_mean

# post            4.15
# pre             7.97



## March 10th Goals
To Do:
1. Go through runs that failed and see why
2. Sort out how depth and quality stats work in ANGSD
3. Run depth, missingness, and total SNP stats to compare two types of alignment

# failed jobs
Slurm Job_id=9276935 Name=angsd_makebeagle Failed, Run time 1-19:46:44, OUT_OF_MEMORY

Slurm Job_id=9276964 Name=angsd_depth Failed, Run time 00:00:12, FAILED, ExitCode 127



# Next steps:
1. Incorporate first low coverage sequence data for post-samples from Sabrina
DONE
2. Split by species
3. Filter by quality >100 
     # justin uses a score of > 20
4. Filter by depth of 4
5. Use 80% missingness as our filtering step
Keep notes on how many snps get kept in each step.
Repeat analyses with depth of 3 and 60% missingness.



sed 's/\SRR2917/SRR/g' darwinfinches_snps_multiallelic_allreads.vcf > darwinfinches_snps_multiallelic_allreads.recode.vcf

#!/bin/bash

#SBATCH --job-name=split_for
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_for.%j

module load bcftools

# for # 18
bcftools view -s SM1156,SM1083,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_snps_multiallelic_allreads.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_snps_multiallelic_allreads.vcf


# par # 19

#!/bin/bash

#SBATCH --job-name=split_par
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_par.%j

module load bcftools

bcftools view -s lamich_PARV1,lamich_PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_snps_multiallelic_allreads.recode.vcf  --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_snps_multiallelic_allreads.vcf

# cra # 14

#!/bin/bash

#SBATCH --job-name=split_cra
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=15gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_cra.%j

module load bcftools

bcftools view -s JP4481,JP5410,JP9655,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_snps_multiallelic_allreads.recode.vcf  --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_snps_multiallelic_allreads.vcf

sbatch split_for.sh 
Submitted batch job 1859116
batch split_cra.sh 
Submitted batch job 1859117
sbatch split_par.sh 
Submitted batch job 1859118




# cra 
# filter by quality > 100

bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_snps_multiallelic_allreads.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality.vcf

#  Filter by depth of 4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4

#  Use 80% missingness as our filtering step

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4_missing8

# generate individual missingness stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4_missing8.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4_missing8

# generate depth stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4_missing8.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra_quality_dp4_missing8


# for

# filter by quality > 100

bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_snps_multiallelic_allreads.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality.vcf

#  Filter by depth of 4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4

#  Use 80% missingness as our filtering step

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4_missing8

# generate individual missingness stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4_missing8.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4_missing8

# generate depth stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4_missing8.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for_quality_dp4_missing8

# par


# filter by quality > 100

bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_snps_multiallelic_allreads.vcf  > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality.vcf

#  Filter by depth of 4

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality.vcf --minDP 4 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4

#  Use 80% missingness as our filtering step

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4.recode.vcf --max-missing 0.8 --recode --recode-INFO-all --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4_missing8

# generate individual missingness stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4_missing8.recode.vcf --missing-indv --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4_missing8

# generate depth stats

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4_missing8.recode.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par_quality_dp4_missing8

sbatch split_cra.sh 
Submitted batch job 9311619


sbatch split_for.sh 
Submitted batch job 9311623

sbatch split_par.sh 
Submitted batch job 9311624



/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/danny_fastqs/sorted_marked_bams/

/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/sorted_marked_bams/


sbatch variantcalling.slurm 
Submitted batch job 1860162


samtools merge lamich_PL15.merged.bam SRR1607532.sorted.bam SRR1607533.sorted.bam 

ref="/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
ID="darwinfinches"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP SM1156.sorted.marked.bam lamich_PARV1.sorted.marked.rehead.bam	 | bcftools call -mv -V indels > "$ID"_whyisitnotworking.vcf

lamich_PARV1.sorted.marked.rehead.bam	   SM1156.sorted.marked.bam
lamich_PARV1.sorted.marked.rehead.bam.bai  SM1156.sorted.marked.bam.bai

ref="/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
bamdir="/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/sorted_marked_bams/"
ID="darwinfinches"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.sorted.marked.bam | bcftools call -mv -V indels > "$ID"_whyisitnotworking.vcf



# trying a different step first 

MergeBamAlignments


java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MergeSamFiles \
INPUT= SRR1607505.sorted.bam I= SRR1607504.sorted.bam \
OUTPUT= sorted_marked_bams/test/lamich_PARV1.merged.marked.bam 

java -jar /opt/ohpc/pub/apps/picard/2.23.4/libs/picard.jar MarkDuplicates \
INPUT= lamich_PARV1.merged.marked.bam \
OUTPUT= lamich_PARV1.merged.sorted.marked.bam \
METRICS_FILE= lamich_PARV1.duplicate.metrics.txt \
MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 &


finch=(lamich_PARV1)

for u in "${finch[@]}"; do
samtools view -H ${u}.merged.sorted.marked.bam > header.${u}.sam
done

for u in "${finch[@]}"; do
samtools reheader header.${u}.sam ${u}.merged.sorted.marked.bam > ${u}.merged.sorted.marked.rehead.bam
done

for u in "${finch[@]}"; do
samtools index ${u}.merged.sorted.marked.rehead.bam
done


ref="/xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa"
bamdir="/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/sorted_marked_bams/"
ID="darwinfinches"
bcftools mpileup -Ou -f "$ref" -a FORMAT/AD,DP,INFO/AD,SP SM1156.sorted.marked.bam lamich_PARV1.merged.sorted.marked.rehead.bam  | bcftools call -mv -V indels > "$ID"_whyisitnotworking.vcf




finch=(lamich_PARV1 lamich_PARV2 lamich_PL15 lamich_PL16 lamich_PL4 lamich_PL7 lamich_PL9)

for u in "${finch[@]}"; do
samtools view -H ${u}.sorted.marked.bam > header.${u}.sam
done

grep 'SM' header.lamich_PARV1.sam

sed -i 's/SM\:SRR1607504/SM\:lamich\_PARV1/g' header.lamich_PARV1.sam
sed -i 's/SM\:SRR1607505/SM\:lamich\_PARV1/g' header.lamich_PARV1.sam

grep 'SM' header.lamich_PARV2.sam
sed -i 's/SM\:SRR1607506/SM\:lamich\_PARV2/g' header.lamich_PARV2.sam
sed -i 's/SM\:SRR1607507/SM\:lamich\_PARV2/g' header.lamich_PARV2.sam

grep 'SM' header.lamich_PL15.sam
sed -i 's/SM\:SRR1607532/SM\:lamich\_PL15/g' header.lamich_PL15.sam
sed -i 's/SM\:SRR1607533/SM\:lamich\_PL15/g' header.lamich_PL15.sam

grep 'SM' header.lamich_PL16.sam
sed -i 's/SM\:SRR1607534/SM\:lamich\_PL16/g' header.lamich_PL16.sam
sed -i 's/SM\:SRR1607535/SM\:lamich\_PL16/g' header.lamich_PL16.sam
sed -i 's/SM\:SRR1607536/SM\:lamich\_PL16/g' header.lamich_PL16.sam

grep 'SM' header.lamich_PL4.sam
sed -i 's/SM\:SRR1607537/SM\:lamich\_PL4/g' header.lamich_PL4.sam
sed -i 's/SM\:SRR1607538/SM\:lamich\_PL4/g' header.lamich_PL4.sam

grep 'SM' header.lamich_PL7.sam
sed -i 's/SM\:SRR1607539/SM\:lamich\_PL7/g' header.lamich_PL7.sam
sed -i 's/SM\:SRR1607540/SM\:lamich\_PL7/g' header.lamich_PL7.sam

grep 'SM' header.lamich_PL9.sam
sed -i 's/SM\:SRR1607541/SM\:lamich\_PL9/g' header.lamich_PL9.sam
sed -i 's/SM\:SRR1607542/SM\:lamich\_PL9/g' header.lamich_PL9.sam

for u in "${finch[@]}"; do
samtools reheader header.${u}.sam ${u}.sorted.marked.bam > ${u}.sorted.marked.rehead.bam
done

for u in "${finch[@]}"; do
samtools index ${u}.sorted.marked.rehead.bam
done


finch=(JP4481 JP5410 JP9655 RHC097 RHC507 SM031 SM032 SM040 SM059 SM079 SM1067 SM1157 SM1200 SM1204 SM1240 SM1266)

for u in "${finch[@]}"; do
samtools view -H ${u}.all.sorted.marked.bam > header.${u}.sam
done

for u in "${finch[@]}"; do
sed -i "s/SM\:${u}\_round1/SM\:${u}/g" header.${u}.sam
done

grep 'SM' header.JP4481.sam
grep 'SM' header.JP5410.sam
grep 'SM' header.JP9655.sam
grep 'SM' header.RHC097.sam
grep 'SM' header.RHC507.sam
grep 'SM' header.SM031.sam
grep 'SM\:' header.SM032.sam
grep 'SM\:' header.SM040.sam
grep 'SM\:' header.SM059.sam
grep 'SM\:' header.SM079.sam
grep 'SM\:' header.SM1067.sam
grep 'SM\:' header.SM1157.sam
grep 'SM\:' header.SM1200.sam
grep 'SM\:' header.SM1204.sam
grep 'SM\:' header.SM1240.sam
grep 'SM\:' header.SM1266.sam




for u in "${finch[@]}"; do
samtools reheader header.${u}.sam ${u}.all.sorted.marked.bam > ${u}.sorted.marked.rehead.bam
done

for u in "${finch[@]}"; do
samtools index ${u}.sorted.marked.rehead.bam
done

# move all modified bams into new folder along with bams of individuals that weren't modified

mv *rehead* finalbams/
mv lamich* finalbams/
mv SRR* finalbams 

finch=(SM1231 SM1270 SM1271 SM1272 SM1273 SM1083 SM1156 SM1237)

for u in "${finch[@]}"; do
mv ${u}* finalbams/
done


/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/sorted_marked_bams/finalbams

sbatch variantcalling.slurm 
Submitted batch job 1861005



# angsd

~/programs/angsd/angsd -beagle /xdisk/mcnew/dannyjackson/finches/angsty/darwinfinches.beagle.gz -doDepth 1 -doCounts 1 -out depth.txt -fai /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa.fai

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -doCounts 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/counts.all -nThreads 8


sbatch angsd_counts.sh 
Submitted batch job 1866455
# says it completed but i can't find the output file anywhere, just a file with the arguments...

sbatch angsd_makebeagle.sh 
Submitted batch job 1866457
# ran out of memory, increased to 1000gb

sbatch angsd_makebeagle.sh 
Submitted batch job 9361847


# March 18th

The counts command in angsd isn't working. I don't understand how or why. I'd like to watch some tutorial videos on ANGSD today and also start writing out a pipeline for the analyses that I'd like to run in ANGSD rather than just continuing to dick around and not get things done.

Play around with PSMC too... compare pre and post estimates of effective pop size
dadi

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -r NC_044571.1-doDepth 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/depth.alls 1 -nThreads 12

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -r NC_044571.1 -doCounts 1 -doDepth 1 -doQsDist 1 -doMajorMinor 2 -doError 1 -out stats.NC_044571 -nThreads 8 # -setMinChunkSize 100

# killed

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -r NC_044571.1 -doCounts 1 -doDepth 1 -doMajorMinor 2 -out stats.NC_044571 -nThreads 8 



NC_044571.1

interactive -a mcnew -n 12 -t 3-00:00 -m 1000gb



Population genomics of the white-beaked dolphin (Lagenorhynchus albirostris): Implications for conservation amid climate-driven range shifts

(-minMapQ 30, -minQ 30, -SNP_pval 1e−6)

#!/bin/bash

#SBATCH --job-name=angsd_makebeagle
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_makebeagle.%j

module load python 

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/finches_initialcheck -nThreads 10 -minMapQ 30 -minQ 30 -doPost 1 -doMajorMinor 1 -doMaf 2 -doGeno 4 -doPlink 2 -SNP_pval 1e−6 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

# tried to do -doPlink 1 and it gave this error:
# Binary plink is not suported yet. use -doPlink 2 
# Must supply -doGeno  to write plink files (consider supplying negative value for suprresing .geno.gz output)
# ran into several more errors requiring that i include flags for doPost, doMajorMinor, and doMaf. used my best interpretation for which flags to use, several gave errors and said they could only work on genotype frequency outputs not on a workflow for calling genotypes, and for those i defaulted to the flags in the example for doPlink in the angsd documents

sbatch angsd_forplink.sh 
Submitted batch job 1887971


# initial inspection of dataset
# Initial variant calling step: use base call and mapping quality filters (-minMapQ 30, -minQ 30, -SNP_pval 1e−6) and write output to PLINK format by specifying the -doPlink flag

# in progress

# inspect this dataset for levels of missing data and distribution of heterozygosity using the --het and --missing functions within PLINK v.1.09 (Purcell et al. 2007). 
# remove samples with missing data ≥ 30% and three samples or above-average heterozygosity suggesting crosscontamination issues
# calculate pairwise relatedness between individuals by combining output from the PLINK --genome function and output from the programme NGSRELATE (Korneliussen and Moltke 2015). For pairs of samples that show a pairwise relatedness coefficient (PI_HAT) above 0.5 corresponding to first-degree relatedness (parent-offspring or full siblings), remove the sample with the lower genotyping rate of each pair

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/darwinfinches -nThreads 10 -doGlf 2 -doMajorMinor 1  -doMaf 2 -SNP_pval 1e-6  -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt


#!/bin/bash

#SBATCH --job-name=angsd_forpopgen
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_forpopgen.%j

module load python 

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/finches_forpopgen -nThreads 10 -minMapQ 30 -minQ 30 -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -minInd 51 -setMinDepthInd 4 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

sbatch angsd_forpopgen.sh 
Submitted batch job 1888052


# -setMaxDepth 1020 is not working

#!/bin/bash

#SBATCH --job-name=angsd_forpopgen2
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_forpopgen2.%j

module load python 

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/finches_forpopgen2 -nThreads 10 -minMapQ 30 -minQ 30 -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doCounts 2 -setMinDepth 204 -setMaxDepth 1020 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

sbatch angsd_forpopgen2.sh 
Submitted batch job 1888054


# dataset for population analyses 
# use genotype likelihood calculation in ANGSD with additional filters on read depth (-setMinDepth 785,-setMaxDepth 3140) corresponding to a minimum depth of coverage of 5X and a maximum depth of coverage of 20X per locus per individual to avoid potential bases arising from sequencing errors following recommendations by O’Leary et al. (2018). 

# Furthermore, we identified variants that were located in an interspersed repeat region using the programme RepeatMasker and excluded those from the variant calling by specifying the remaining sites using the -sites flag in ANGSD. Further filtering of the multilocus genotypes was conducted in PLINK using a minor allele count of 2 to remove variants generated through uncertainties in base calling during sequencing. We examined the patterns of linkage disequilibrium decay in our data and observed a relatively steep decline in linkage disequilibrium in the initial portion of your linkage disequilibrium decay graph drawn by the programme NGSLD (Fox et al. 2019). This suggests stronger linkage patterns among nearby SNPs and therefore, we used the --indep function in PLINK to prune loci affected by linkage disequilibrium with a window size of 50 kb, a step size of 5 and a variant inflation factor of 2. The final dataset comprised 1092 Single
Nucleotide Polymorphisms (SNPs) for all downstream population genetic
analyses.





# depth stats on all reads

vcftools --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic.vcf --depth --out /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic




sed 's/\SRR2917/SRR/g' darwinfinches_all_snps_multiallelic.vcf > darwinfinches_all_snps_multiallelic.vcf.recode.vcf

# for # 18

#!/bin/bash

#SBATCH --job-name=split_for
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_for.%j

module load bcftools

bcftools view -s SM1156,SM1083,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298  /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic.vcf.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/for.vcf


# par # 19

#!/bin/bash

#SBATCH --job-name=split_par
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_par.%j

module load bcftools

bcftools view -s lamich_PARV1,lamich_PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic.vcf.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/par.vcf

# cra # 14

#!/bin/bash

#SBATCH --job-name=split_cra
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.split_cra.%j

module load bcftools

bcftools view -s JP4481,JP5410,JP9655,lamich_PL15,lamich_PL16,lamich_PL4,lamich_PL7,lamich_PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic.vcf.recode.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/allreads/cra.vcf

# monday march 18 2024 9:58pm
(elgato) [dannyjackson@cpu37 slurmscripts]$ sbatch split_for.sh 
Submitted batch job 1872414
(elgato) [dannyjackson@cpu37 slurmscripts]$ sbatch split_cra.sh 
Submitted batch job 1872415
(elgato) [dannyjackson@cpu37 slurmscripts]$ sbatch split_par.sh 
Submitted batch job 1872416



# angsd 
# pca
https://github.com/Rosemeis/pcangsd

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt -r NC_044571.1 -doCounts 1 -doDepth 1 -doQsDist 1 -doMajorMinor 2 -doError 1 -out stats.NC_044571 -nThreads 8


#!/bin/bash

#SBATCH --job-name=pcangsd
#SBATCH --ntasks=2
#SBATCH --nodes=1             
#SBATCH --time=12:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=500gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pcangsd.%j

pcangsd -b /xdisk/mcnew/dannyjackson/finches/angsty/darwinfinches.beagle.gz -o pca_output -t 8

sbatch pcangsd_pca.sh 
Submitted batch job 1872407


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


# compare to bcftools call depth
python
import numpy as np
import pandas as pd 

df = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/vcfs/allreads/darwinfinches_all_snps_multiallelic.idepth', sep='\t')


df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = df.join(df_names.set_index('sample'), on='INDV')

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')

                      N_SITES  MEAN_DEPTH
species treatment                        
CRA     post       5.795318
        pre        10.330662
FOR     post       4.795200
        pre        6.820301
PAR     post       5.110453
        pre        9.119531

module load htslib
bgzip darwinfinches_all_snps_multiallelic.vcf
bcftools view -r NC_044571.1 darwinfinches_all_snps_multiallelic.vcf > NC_044571.vcf


# bcftools view -r NC_044571.1,NC_044572.1,NC_044573.1,NC_044574.1,NC_044575.1,NC_044576.1,NC_044577.1,NC_044578.1,NC_044579.1,NC_044580.1,NC_044581.1,NC_044582.1,NC_044583.1,NC_044584.1,NC_044585.1,NC_044586.1,NC_044587.1,NC_044588.1,NC_044589.1,NC_044590.1,NC_044591.1,NC_044592.1,NC_044593.1,NC_044594.1,NC_044595.1,NC_044596.1,NC_044597.1,NC_044598.1,NC_044599.1,NC_044600.1,NC_044601.1 /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par.autosomes.vcf



# March 19th plan for filtering

Main analyses:
Don't drop any individuals 
Filter stringently
Use ANGSD

Sanity check analyses:
Drop individuals to maintain similar depth scores
Filter lightly
Use bcftools



# pcangsd 
module load python
pcangsd -b /xdisk/mcnew/dannyjackson/finches/angsty/darwinfinches.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/angsty/PCA/pca_output -t 8


python3 $PCANGSD -beagle $BASEDIR/results/MME_ANGSD_PCA_LDpruned.beagle.gz -o $BASEDIR/results/PCAngsd_LDpruned_covmat

sbatch pcangsd_pca.sh 
Submitted batch job 9382749

pop_code <- read.csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')





library(ggfortify)

name <- "pca_output.cov"
m <- read.table(name)
pop_code <- read.csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

colnames(m) <- pop_code$sample
rownames(m) <- pop_code$sample

pca_res <- prcomp(m, scale. = TRUE)


pdf(file = "pca_species.pdf", width = 10, height = 10, useDingbats=FALSE)
autoplot(pca_res, data = pop_code, color = 'species')
dev.off()

pdf(file = "pca_indiv.pdf", width = 10, height = 10, useDingbats=FALSE)
autoplot(pca_res, data = pop_code, color = 'species', label=TRUE)
dev.off()


# Thursday March 21st
Main analyses:
Don't drop any individuals 
Filter stringently
Use ANGSD

Sanity check analyses:
Drop individuals to maintain similar depth scores
Filter lightly
Use bcftools

Goals for today:
1. Write ANGSD code with proper filters
2. Filter vcf from bcftools
  I don't think I actually need to do this... 
3. 



# Friday March 22nd
Reading up on batch effects
Use Trimmomatic SLIDING WINDOW to drop base calls with a low quality at the end of reads
This approach is implemented in fastp as the cut_right op- tion, and in trimmomatic as the SLIDINGWINDOW option. Because a drop in base quality is often not immediately followed by a poly-G tail, sliding-window base quality trimming may result in greater data loss than necessary, but we found it to be much more effective at removing poly-G tails than targeted poly-G trimming with existing tools (Figure 2a). Indeed, after applying this method (with window size of 4 and average base quality threshold of 20), G bases are no longer enriched at the end of reads in our samples sequenced in the NextSeq-150PE batch (Figure 2b). 

# Add in Trimmomatic trimming step 
# I followed sabrina's code which just used adapterremoval-2.3.1
# AdapterRemoval --file1 ${u}_1.fastq.gz --file2 ${u}_2.fastq.gz --trimns --trimqualities --minquality 20 --minlength 25 --collapse --threads 8  --basename ${u}

# Use a more stringent base quality filter of 33 (rather than 20)




# mapping filters can differentially affect reads of different lengths, resulting in false outlier loci between sequencing types
# Based on this pattern, a simple mitigation strategy is to locate the sites that have a high proportion of low-mapping-score reads mapping to them (e.g., >10%) in a batch of data with single-end reads and/or shorter reads and exclude them from further analyses.


# DNA degradation can lead to deamination of cytosines, and an increased frequency of C-T transversions, also resulting in higher estimates of changes in nucleotide diversity between popualtions of different time periods.
# the change in diversity estimates after excluding all C-to-T and G-to-A transitions (e.g., the -noTrans 1 option in angsd. Ignoring a subset of variant types certainly results in decreases in diversity indices in all samples, but if some samples are more strongly impacted, it means that DNA deg- radation levels are uneven among samples (Figure 5d). In this case, the diversity estimates excluding transitions will be more compara- ble between batches and less biased in a relative sense (Figure 5a “excluding transitions”). We also observed that after reference bias is corrected for, there is no notable batch effects in our PCA results when transitions are included (Figure 5b “including transitions”), and excluding transitions does not make a significant difference (Figure 5b “excluding transitions”), suggesting that the deamination of cytosines is not a major cause of batch effects in PCA for our data (but it could be in other data sets where the samples have suffered greater postmortem damage).

# adjusting for differences in read depth for PCA analyses
# The only differ- ence between the two batches of simulated data is their sequenc- ing depth (either 0.125 or 4×, see Supporting Information for details about the simulations). At low sample size (five or 10 per population), PCAs generated from pcangsd-0.98 (Meisner & Albrechtsen, 2018) and the -doCov 1 option in angsd tend to group samples with the same read depth together along one of the top PC axes, creating false patterns of clustering (Figure 6). In comparison, the PCoA gen- erated from the -doIBS 2 option in angsd is less prone to such biases (Figure 6). We observed a similar pattern in our empirical data, where the PCA generated from the -doCov 1 option in angsd does not show obvious signs of batch effects when other causes of batch effects are controlled for, despite the difference in sequencing depth be- tween the two batches (Figure 1b “after” and 5b). In contrast, PCA generated from pcangsd still has individuals from different batches clustering separately (Figure S6).

# However, if samples are not randomly assigned and if true bio- logical signals may be confounded with batch effects, it may no lon- ger be possible to determine the presence/absence of batch effects. In such cases, we would recommend researchers to take a subset of data from each batch, and perform some of the tests that we have mentioned in this paper (comparing heterozygosity estimates be- fore and after applying a stringent base quality filter, calculating the frequencies of different base substitutions in private alleles in each batch of data, etc.) as a means to determine the presence/absence of batch effects.

# Disparities in sequencing depth is unlikely to become an issue if depth is higher than 20× in all batches. Otherwise, genotype calling in the batch with lower coverage (even at medium coverage, e.g., 5–20×) is likely to be more inaccurate, and may therefore cause batch effects (Warmuth & Ellegren, 2019). In these cases, genotype- likelihood-based inference may be preferable to genotype calling. 





# Run tests of biases

Testing for and mitigating various types of batch effects
1. Test for Batch Effects of Depth
Question: 
  Does depth of coverage affect population statistics? 

How to test for presence of batch effects:
  Downsample the batch of data with higher coverage and compare the results generated from before and after downsampling

Recommended mitigation of batch effects: 
  PCA generated from the -doCov 1 option in angsd does not show obvious signs of batch effects when other causes of batch effects are controlled for, despite the difference in sequencing depth between the two batches (Figure 1b “after” and 5b). In contrast, PCA generated from pcangsd still has individuals from different batches clustering separately (Figure S6).



2. Test for Batch Effects of Age / degradation

How to test for presence of batch effects:
  1. Compare the frequencies of different types of base substitutions among the private alleles in each batch (Figure 5c)
  2. Compare the drop in diversity estimates (e.g., individual heterozygosity) after excluding all transitions between different batches of data (Figure 5d)

Recommended mitigation of batch effects: 
  1. Exclude transitions from certain analyses (Figure 5a “excluding transitions”; Figure S3)
  2. Recalibrate base quality scores for degradated DNA (e.g., mapdamage) (Jónsson et al., 2013)
  3. Use genotype likelihood models that take postmortem damage into account (e.g., atlas) (Link et al., 2017)



3. Test for Batch Effects of Sequencing platform 
Question 1: 
  Do our samples differ in their sequencing platform?
How to test for presence of batch effects: 
  Look up sequencing platform used in other sequencing

Question 2: 
  Do sequencing platforms affect population statistics?

How to test for presence of batch effects:
  Examine the base composition at each read position in raw fastq files (e.g., with fastqc)

Recommended mitigation of batch effects: 
  Trim off ends of reads with low base quality within sliding windows (e.g., the cut_right option in fastp, or the SLIDINGWINDOW option in trimmomatic)



4. Test for Batch Effects of Using Duplicate Sequencing of Individuals
Question: 
  Does combining sequencing runs of the same individual to increase depth lead to different levels of miscalibration in base quality scores?

# I'm treating this as the same question as "Do separate sequencing runs lead to different levels of miscalibration in base quality scores?"

How to test for presence of batch effects:
  Compare diversity estimates(e.g.,individual heterozygosity) using a relaxed vs. stringent base quality threshold within each batch

Recommended mitigation of batch effects: 
  1. Use a more stringent base quality threshold in all batches (Figure 1a “after” and 3; Figure S3)
  2. Use base quality score recalibration (Figure S4) if there is a comprehensive variant database (e.g., the soapSNP model in angsd, gatk) (Korneliussen et al., 2014; McKenna et al., 2010) or if control sequences (e.g., PhiX) are included in the sequencing run (Ni & Stoneking, 2016; Orr, 2020; Zook et al., 2012)



5. Test for Batch Effects of Sample size
Question:
  Does uneven sample size affect population statistics?
  
How to test for presence of batch effects:
  Not outlined in the Nina papers.
  I guess just compare pop stats using all samples and again with even subsets of samples. Could randomize subsets of samples and create a "null distribution" even to compare to the uneven sample size.

Recommended mitigation of batch effects: 
  Idk bro maybe downsample if it isn't also biased








# Start writing code for these different runs
# What filtering parameters should I have as my default?
Platform
Dups
Depth
Age
Sample Size


3. Test for Batch Effects of Sequencing platform 
Question 1: 
  Do our samples differ in their sequencing platform?
How to test for presence of batch effects: 
  Look up sequencing platform used in other sequencing
Unfortunately yes they do differ

Question 2: 
  Do sequencing platforms affect population statistics?

How to test for presence of batch effects:
  Examine the base composition at each read position in raw fastq files (e.g., with fastqc)

Recommended mitigation of batch effects: 
  Trim off ends of reads with low base quality within sliding windows (e.g., the cut_right option in fastp, or the SLIDINGWINDOW option in trimmomatic)



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





















# Monday March 25th
# Things that finished running over the weekend
angsd_forpopgen 1888052
angsd_makebeagle 1887971
angsd_forpopgen2 1888141 # really fastqc_rawfastas.sh, not forpopgen2
angsd_forpopgen2 1888054

# angsd_forpopgen2 1888141 # really fastqc_rawfastas.sh, not forpopgen2
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

# Illumina universal adapter is present in many samples
# GC content is flagged in almost if not all samples
# Trimmomatic then repeat

# base filtering steps:
NextSeq-150PE: PE -phred33 'ILLUMINACLIP:'$ADAPTERS':2:30:10:1:true', HiSeq-125SE: SE -phred33 'ILLUMINACLIP:'$ADAPTERS':2:30:10'

used fastp-0.19.7 to trim poly-G tails with the NextSeq-150PE batch of data only (--trim_poly_g -Q -L -A), with the default setting on minimum poly-G length threshold (--poly_g_min_len 10

# output fastas to here:
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive
# TruSeq3 adapters are used in HiSeq and novaseq: TruSeq3-PE.fa TruSeq3-PE-2.fa 
/home/u15/dannyjackson/programs/Trimmomatic/adapters/TruSeq3-PE-2.fa 
/home/u15/dannyjackson/programs/Trimmomatic/adapters/TruSeq3-PE.fa 

# shell script for pre-trim QC, trimming, and post-trim QC
# pre and post trim QC files go to separate directories
cat /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames.txt > /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt

ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*round1.sorted.bam | awk 'BEGIN {FS = "/"} {print $7}' | awk 'BEGIN {FS = "."} {print $1}' >> /xdisk/mcnew/dannyjackson/finches/reference_lists/sortedbam_samplenames.txt

ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1*gz | awk 'BEGIN {FS = "/"} {print $7}' | awk 'BEGIN {FS = "_"} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt

ls /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/*_1*gz | grep round1 | awk 'BEGIN {FS = "/"} {print $7}' | awk 'BEGIN {FS = "_"} {print $1,".round1"}' | sed 's/ //g' > /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

awk 'BEGIN {FS = ",";OFS = ","} $3="NC_0"$3' test.csv > test2.csv


#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=60:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=600gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=trimming.%j

cd ~/programs/Trimmomatic

echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt

while read -r ID;
do
  java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 12 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_1.fastq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_2.fastq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trimmed_fastas/"$ID"_trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:1:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt


done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt

sbatch ~/programs/slurmscripts/trim.sh 
Submitted batch job 1889099

Exception in thread "main" java.io.FileNotFoundException: /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SM1266_round1_1.fastq.gz (No such file or directory)

java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 12 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SM1266_1.round1.fq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SM1266_2.round1.fq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trimmed_fastas/SM1266_1_trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:1:30:10 \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt


# batch effect trimming steps:
used the sliding window quality trimming functionality in fastp-0.19.7 (--cut_right --trim_poly_g -L -A) to further eliminate poly-G tails in the adapter trimmed fastq files in the NextSeq-150PE batch. Default window length (--cut_right_window_size 4) and mean base quality threshold (--cut_right_mean_quality 20) 
Could instead use the sliding window function in trimmmomatic



# angsd_makebeagle 1887971
finches_initialcheck.tped

module load plink

/xdisk/mcnew/dannyjackson/finches/angsty/finches_initialcheck.tped

# inspect this dataset for levels of missing data and distribution of heterozygosity using the --het and --missing functions within PLINK v.1.09 (Purcell et al. 2007). 
# remove samples with missing data ≥ 30% and three samples or above-average heterozygosity suggesting crosscontamination issues
# calculate pairwise relatedness between individuals by combining output from the PLINK --genome function and output from the programme NGSRELATE (Korneliussen and Moltke 2015). For pairs of samples that show a pairwise relatedness coefficient (PI_HAT) above 0.5 corresponding to first-degree relatedness (parent-offspring or full siblings), remove the sample with the lower genotyping rate of each pair

#!/bin/bash

#SBATCH --job-name=plink_missingandhet
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=60:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

module load plink

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/developing_pipeline

plink --tped /xdisk/mcnew/dannyjackson/finches/angsty/finches_initialcheck.tped --tfam /xdisk/mcnew/dannyjackson/finches/angsty/finches_initialcheck.tfam --allow-extra-chr --missing --freq --het

sbatch ~/programs/slurmscripts/plink_missingandhet.sh 
Submitted batch job 9424697


# # Tuesday March 26th
# Things that finished running over the weekend that I haven't addressed
angsd_forpopgen 1888052
angsd_forpopgen2 1888054

# Things that finished running overnight

# evaluate trimmmomatic output
# run another trimming step
# Look at both forpopgen runs and continue building out their pipelines

Trimming is done, so run fastqc and reevaluate stats post-trim

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

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/platform/trimmed_fastqcs

fastqc -t 12 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trimmed_fastas/*.fq.gz

sbatch ~/programs/slurmscripts/faimmed.sh 
Submitted batch job 9426128


# evaluating the forpopgen and forpopgen2 outputs


#!/bin/bash

#SBATCH --job-name=angsd_forpopgen2
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_forpopgen2.%j

module load python 

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/finches_forpopgen2 -nThreads 10 -minMapQ 30 -minQ 30 -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doCounts 2 -setMinDepth 204 -setMaxDepth 1020 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

sbatch angsd_forpopgen2.sh 
Submitted batch job 1888054


#!/bin/bash

#SBATCH --job-name=angsd_forpopgen
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=angsd_forpopgen.%j

module load python 

~/programs/angsd/angsd -GL 1 -out /xdisk/mcnew/dannyjackson/finches/angsty/finches_forpopgen -nThreads 10 -minMapQ 30 -minQ 30 -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -minInd 51 -setMinDepthInd 4 -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/bamlist.txt

# initialcheck plink output
Relatively similar and low imiss numbers, but really variable (and likely biased by sequencing) heterozygosity scores. I'm stuck a little on this step because the trimming should address a lot of these errors. The pipeline is there and I can run it when things are a little more in line. Wait for the fastqcs to wrap. Consider trimming for polyg tails before using slidingwindow in trimmomatic... this is exhausting the number of potential pathways that I'm evaluating.

scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/finches/bias_testing/platform/trimmed_fastqcs/ .


(--trim_poly_g -Q -L -A), with the default setting on minimum poly-G length threshold (--poly_g_min_len 10). 


cmake --build build -prefix "/home/u15/dannyjackson/programs/libdeflate-1.19"


cmake --install build -prefix "/home/u15/dannyjackson/programs/libdeflate-1.19"

cmake --install . --prefix "/home/u15/dannyjackson/programs/libdeflate-1.19"


make install

export PATH=/home/u15/dannyjackson/programs/libdeflate-1.19:$PATH
export PATH=/home/u15/dannyjackson/programs/isa-l:$PATH



#!/bin/bash

#SBATCH --job-name=trimming.%j
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

cd ~/programs/Trimmomatic

echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt

while read -r ID;
do

  echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/adapter_trim_log.txt

  java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 24 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_1.fastq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_2.fastq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10:1 >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/trim_log.txt


done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt

while read -r ID;
do

echo "Beginning polyg trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_1P.fq.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_2P.fq.gz -O /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz 

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_samplenames.txt

sbatch ~/programs/slurmscripts/trim_adapters_polyg.sh 
Submitted batch job 9440926

sbatch trim_adapters_polyg.sh 
Submitted batch job 9447584


# error occured in the polyg script, rerunning with just that step. see Friday March 29th notes

# Thursday March 28th
I'm now running the proper script for trimming adapters followed by polyg trimming. I can spend my time writing the script for trimming with sliding window filters after poly g trimming, and for generating fastq files from each step of filtering. I'd also like a more straightforward output to read than the htmls but maybe that isn't worth my time.


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

sbatch trim2.sh 
Submitted batch job 9477536

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
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_round1.polygtrimmed.fastq.gz  /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_round1.polygtrimmed.fastq.gz  \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/"$ID"_round1.trimmed2.fq.gz \
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim2_log.txt


done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

sbatch trim2_round1.sh 
Submitted batch job 9477560

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

# output warn and fail flags from fastqc data
# unpack zip files and navigate to directory above all zip files
grep 'warn' */fastqc_data.txt
grep 'fail' */fastqc_data.txt

# angsd notes

# intial pass
snp filter, maf filter (jesse prefers not), depth filter, include minor allele frequency > output global snp list

If it keeps breaking on memory, might be best to separate by chromosome.

Can create either a beagle or a vcf, both seem to hold freq data. Can also directly compute a site frequency spectrum that incorporates uncertainty 

Creates a list of all the sites so that downstream we only consider these sites

Nicolas (lcwgs papers) has a dissertation on historical and contemporary atlantic cod

Nina has a paper in Science (diff stat method normalizes differences in allele frequencies)

Filter at species level not with all 
Jesse filters less stringently (hesitant to filter on minor allele freq, missingness; does filter on depth)



GATK indel realignment: look into this, it's a rabbit hole. There may have been a shift between GATK 2 vs 3... it may not be simple.
gatk/4.2.5.0
After GATK 3.6, no need to do indel realignment because HAPLOTYPECALLER does it on the fly.





#!/bin/bash

#SBATCH --job-name=trimming
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=240:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=1000gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

cd ~/programs/Trimmomatic

echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/trim_log.txt

while read -r ID;
do

  echo "Beginning trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/adapter_trim_log.txt

  java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 24 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_1.fastq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_2.fastq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10:1 >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/trim_log.txt


done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

while read -r ID;
do

echo "Beginning polyg trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_trimmed_1P.fq.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_polygtrimmed.fastq.gz -I /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_trimmed_2P.fq.gz -O /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_polygtrimmed.fastq.gz 

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt


#!/bin/bash

#SBATCH --job-name=polyg_trimming
#SBATCH --ntasks=24
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=500gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.%j

while read -r ID;
do

echo "Beginning polyg trimming for "$ID>>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/polyg_trim_log.txt

  ~/programs/fastp --trim_poly_g -Q -L -A -w 24 --poly_g_min_len 10 -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"round1.trimmed_1P.fq.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_1_round1.polygtrimmed.fastq.gz -I /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/"$ID"_round1.trimmed_2P.fq.gz -O /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/polygtrimmed_fastas/"$ID"_2_round1.polygtrimmed.fastq.gz 

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

sbatch polyg_trimming_round1.sh 
Submitted batch job 9456484

sed -i 's/\.round1//g' /xdisk/mcnew/dannyjackson/finches/reference_lists/fasta_round1_samplenames.txt

sbatch trim_adapters_polyg_round1.sh 
Submitted batch job 9456470

/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SM1266_1.round1.fastq.gz
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/"$ID"_1.round1.fastq.gz

# Friday March 29th 

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




# Monday April 1st
Running fastqc on adapter trimmed files

/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/

# redo anything that's less than 1G with an interactive session
-rw-r--r--. 1 dannyjackson mcnew 711M Mar 29 14:58 SRR1607535_trimmed_2P.fq.gz
-rw-r--r--. 1 dannyjackson mcnew 918M Mar 31 13:14 JP9655_round1.trimmed_2P.fq.gz

#!/bin/bash

#SBATCH --job-name=trim_test
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=60:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=500gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=trim_test.%j

cd ~/programs/Trimmomatic

java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 6 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SRR1607535_1.fastq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/SRR1607535_2.fastq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/SRR1607535_test_round1.trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10:1 >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/trim_log.txt

java -jar ~/programs/Trimmomatic/dist/jar/trimmomatic-0.40-rc1.jar PE -threads 6 \
/xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/JP9655_1.round1.fq.gz /xdisk/mcnew/dannyjackson/finch_wgs/danny_fastqs/JP9655_2.round1.fq.gz \
-baseout /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/JP9655_test_round1.trimmed.fq.gz \
ILLUMINACLIP:adapters/TruSeq3-PE-2.fa:2:30:10:1 >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/adaptertrimmed_fastas/trim_log.txt

sbatch trim_test.sh 
Submitted batch job 1895737

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


#!/bin/bash

#SBATCH --job-name=fastqc_trimmed2
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=60:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.#SBATCH --job-name=fastqc_trimmed.%j

module load fastqc/0.11.9

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/

fastqc -t 12 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/finaltrim_fastas/*.fq.gz

sbatch fastqc_trim2.sh 
Submitted batch job 9480270

grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_adaptertrimmed_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_adaptertrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' failings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_polygtrimmed_fastas.txt

awk 'BEGIN {FS = ">>"} {print $2}' warnings.txt | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_polygtrimmed_fastas.txt

grep '_2P_' failings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_fail_finaltrim_fastas.txt

grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_finaltrim_fastas.txt



# fastqc report summaries, includes all fastas
stats_warn_adaptertrimmed_fastas.txt (726 total)
     56 Overrepresented sequences
      5 Per base N content
    108 Per base sequence content
    182 Per sequence GC content
     17 Per sequence quality scores
     43 Per tile sequence quality
     18 Sequence Duplication Levels
    297 Sequence Length Distribution
cat stats_fail_adaptertrimmed_fastas.txt (503 total)
     62 Overrepresented sequences
     91 Per base sequence content
     85 Per base sequence quality
    114 Per sequence GC content
     35 Per sequence quality scores
    107 Per tile sequence quality
      9 Sequence Duplication Levels

stats_warn_polygtrimmed_fastas.txt (389 total)
      6 Overrepresented sequences
     59 Per base sequence content
    133 Per sequence GC content
     41 Per tile sequence quality
    150 Sequence Length Distribution
stats_fail_polygtrimmed_fastas.txt (17 total)
      2 Per base sequence quality
     13 Per sequence GC content
      2 Per tile sequence quality

cat stats_warn_finaltrim_fastas.txt (685 total)
    100 Per base sequence content
    240 Per sequence GC content
     45 Per tile sequence quality
    300 Sequence Length Distribution
cat stats_fail_finaltrim_fastas.txt (61 total)
      1 Per base sequence content
     57 Per sequence GC content
      3 Per tile sequence quality

# fastqc report summaries, only 2P fastas
stats_warn_adaptertrimmed_fastas.txt 
      6 Overrepresented sequences
     25 Per base sequence content
     66 Per sequence GC content
     20 Per tile sequence quality
     75 Sequence Length Distribution
stats_fail_adaptertrimmed_fastas.txt 
      1 Per base sequence quality
      6 Per sequence GC content
      1 Per tile sequence quality

stats_warn_polygtrimmed_fastas.txt (389 total)
      6 Overrepresented sequences
     59 Per base sequence content
    133 Per sequence GC content
     41 Per tile sequence quality
    150 Sequence Length Distribution
stats_fail_polygtrimmed_fastas.txt (17 total)
      2 Per base sequence quality
     13 Per sequence GC content
      2 Per tile sequence quality

stats_warn_finaltrim_fastas.txt 
     26 Per base sequence content
     67 Per sequence GC content
     75 Sequence Length Distribution
stats_fail_finaltrim_fastas.txt 
      7 Per sequence GC content
      1 Per tile sequence quality




We didn’t improve on three metrics:
1.	Per base sequence content (25 to 26)
2.	Per sequence GC content (72 to 74)
3.	Sequence Length Distribution (75 to 75)

grep '_2P_' warnings.txt | grep 'Per base sequence content' | awk 'BEGIN {FS = "."} {print $1}' | sort | uniq 
grep '_2P_' failings.txt | grep 'Per base sequence content' | awk 'BEGIN {FS = "_"} {print $1}' | sort | uniq 

grep '_2P_' failings.txt | grep 'Per sequence GC content' | awk 'BEGIN {FS = "."} {print $1}' | sort | uniq 


grep '_2P_' warnings.txt | awk 'BEGIN {FS = ">>"} {print $2}' | awk 'BEGIN {FS = "\t"} {print $1}' | sort | uniq -c > stats_warn_finaltrim_fastas.txt



Next Steps:
1. Clip overlapping read pairs
2. Indel realignment (see if this is necessary)
3. Mark duplicates
4. Align
5. Compute likelihoods
6. Call genotypes?

deduplication, overlap clipping, and indel-realignment 

Mapped reads to the reference genome 
  bowtie2-2.3.5.1 (-q --phred33 -- very-sensitive -I 0 -X 1500 –fr for NextSeq-150PE and -q --phred33 -- very-sensitive for HiSeq-125SE). 
  Used samtools-1.11 to convert the resulting sam files to bam format and sorted them (view -buS and sort),

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



 Picard tools-2.9.0 to remove duplicated reads (MarkDuplicates VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true), 
 
sbatch markduplicates.sh 
Submitted batch job 9515765

 
 and BamUtil-1.0.14 to clip overlapping read pairs (clipOverlap) with the NextSeq-150PE batch of data only. 

ls /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/

ls *all*bam | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.all.txt

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

/home/u15/dannyjackson/programs/bamUtil-master/bin/bam clipOverlap --in /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam --out /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam --stats --params


echo "done " $@ >>/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/clippingstats.txt 
}


export -f clipping 

parallel -j 12 clipping :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch clipoverlap_parallel.sh 
Submitted batch job 1934390

# redoing it with the .all. bamfiles
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


 Lastly, we performed indel realignment with GATK-3.7 (-T RealignerTargetCreator followed by -T IndelRealigner --consensusDeterminationModel USE_READ, with default options). 
 
apptainer pull docker://broadinstitute/gatk3:3.6-0

 apptainer pull docker://broadinstitute/gatk3:3.7-0
 
  apptainer exec ~/programs/gatk3_3.7-0.sif  java -jar /usr/GenomeAnalysisTK.jar --help
 apptainer exec ~/programs/gatk3_3.7-0.sif  java -jar /usr/GenomeAnalysisTK.jar -T RealignerTargetCreator 


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

 Lastly, we counted the number of bases with mapping quality higher than 20 in the indel-realigned bam files using Samtools, and calculated per-sample sequencing depth (Table 1, Figure S1). 

samtools depth -q 30 JP4481.realigned.bam

module load samtools
while read -r finch;
do 
  samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

ls *bam | awk 'BEGIN {FS = "."} {print $1}' | wc -l > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_indelrealigned.txt

# test for effects of duplicate sequencing effort
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_round1.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_all.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266.realigned.bam
/xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_round1.realigned.bam


#!/bin/bash

#SBATCH --job-name=genolike
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolike.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist -r 	NC_044571.1

sbatch genolike.chr1.sh 
Submitted batch job 9754424

#!/bin/bash

#SBATCH --job-name=genolike
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolike.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist

slurmscripts]$ sbatch genolike.sh 
Submitted batch job 9754425


/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.chr1.beagle.gz

python pcangsd.py -h



# trying to install pcangsd with cond
git clone https://github.com/Rosemeis/pcangsd.git
cd pcangsd
conda env create -f environment.yml

conda activate pcangsd
pip3 install .
python setup.py build_ext --inplace


export PYTHONPATH="${PYTHONPATH}:/home/u15/dannyjackson/.local/lib/python3.12/site-packages"

/home/u15/dannyjackson/.local/lib/python3.12/site-packages

conda run -n pcangsd python --version




#!/bin/bash

#SBATCH --job-name=pcangsd
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=50:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.pcangsd.%j

conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/output -t 12



# skipping to heterozygosity while waiting for hpc support to help me out with pcangsd

~/programs/angsd/angsd - -i my.bam -anc ref.fa -dosaf 1 -fold 1



#!/bin/bash

#SBATCH --job-name=hetero
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.hetero.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/heterozygosity

~/programs/angsd/angsd -dosaf 1 -out chr1.sfs -nThreads 10 -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist -minQ 30 -minMapQ 30 -setMinDepthInd 4 -r NC_044571.1 -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -gl 1

sbatch angsd_heterozygosity_dups_chr1.sh 
Submitted batch job 1950572



#!/bin/bash

#SBATCH --job-name=hetero
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.hetero.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/heterozygosity

conda activate pcangsd

#followed by the actual estimation
~/programs/angsd/misc/realSFS chr1.sfs.saf.idx > chr1.est.ml


angsd_het_compuation.sh 
Submitted batch job 9792974

angsd_het_compuation.sh 
Submitted batch job 9795857


#in R
a<-scan("chr1.est.ml")
a[1]/sum(a)



C <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


# Plot PCA plot
e <- eigen(C)
pdf(file = "pca_chr1.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

dev.off()




# test for effects of duplicate sequencing effort
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP4481_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP5410_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/JP9655_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1067_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1157_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1200_round1.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_all.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266.realigned.bam
CRA /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_round1.realigned.bam

#!/bin/bash

#SBATCH --job-name=genolikecra
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikecra.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.cra.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.cra -r 	NC_044571.1

sbatch genolike.cra.chr1.sh 
Submitted batch job 9767666

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.cra.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/dups/output -t 12 --iter 10000

awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.cra | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/dups/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/dups/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_cra.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()

# Plot PCA plot
e <- eigen(C)
plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd")

# Obtain p-values from PC-based selection scan
p <- pchisq(D, 1, lower.tail=FALSE)


PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC097_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/RHC507_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM031_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM032_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM040_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM079_round1.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_all.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240.realigned.bam
PAR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_round1.realigned.bam


#!/bin/bash

#SBATCH --job-name=genolikepar
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikepar.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.par.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.par -r 	NC_044571.1

sbatch genolike.par.chr1.sh 
Submitted batch job 9767669

conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.par.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/dups/output -t 12 --iter 10000



awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.par | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/dups/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/dups/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_par.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()


FOR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204_all.realigned.bam
FOR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204.realigned.bam
FOR /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1204_round1.realigned.bam

#!/bin/bash

#SBATCH --job-name=genolikefor
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikefor.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.for.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.for -r 	NC_044571.1

sbatch genolike.for.chr1.sh 
Submitted batch job 9767675



conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.for.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/dups/output -t 12 --iter 10000



awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.for | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/dups/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/dups/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_for.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()


ls /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/*bam > bam.filelist

#!/bin/bash

#SBATCH --job-name=genolike
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolike.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist -r NC_044571.1

sbatch genolike.chr1.sh 
Submitted batch job 9767690

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/
cd all

conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/all/output -t 12 --iter 10000



awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/all/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/all/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_all.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()


# when these are done, do a pca on cra, for, and par separately and one on all bam files 

# this is the list of all bams that are involved in duplicates -- in pca of all bams, color these together and label all as individuals. Should actually additionally run all cra and all for and all par not just dups

JP4481_all
JP4481
JP4481_round1
JP5410_all
JP5410
JP5410_round1
JP9655_all
JP9655
JP9655_round1
RHC097_all
RHC097
RHC097_round1
RHC507_all
RHC507
RHC507_round1
SM031_all
SM031
SM031_round1
SM032_all
SM032
SM032_round1
SM040_all
SM040
SM040_round1
SM079_all
SM079
SM079_round1
SM1067_all
SM1067
SM1067_round1
SM1157_all
SM1157
SM1157_round1
SM1200_all
SM1200
SM1200_round1
SM1204_all
SM1204
SM1204_round1
SM1240_all
SM1240
SM1240_round1
SM1266_all
SM1266
SM1266_round1

# SM059 is missing "all" file

#!/bin/bash

#SBATCH --job-name=genolikecraall
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikecraall.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.cra.all.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.cra -r 	NC_044571.1

sbatch genolike.cra.all.chr1.sh 
Submitted batch job 9767741

conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.cra.all.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/all/output -t 12




awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.cra | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/all/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/cra/all/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_cra_all.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()



#!/bin/bash

#SBATCH --job-name=genolikeforall
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikeforall.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.for.all.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.for -r 	NC_044571.1

sbatch genolike.for.all.chr1.sh 
Submitted batch job 9767740


conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.for.all.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/all/output -t 12 --iter 10000


awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.for | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/all/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/for/all/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_for_all.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()


#!/bin/bash

#SBATCH --job-name=genolikeparall
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.genolikeparall.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/

~/programs/angsd/angsd -GL 2 -out genolike.par.all.chr1 -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.par -r 	NC_044571.1

sbatch genolike.par.all.chr1.sh 
Submitted batch job 9767739


conda activate pcangsd

pcangsd -b /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/beagles/genolike.par.all.chr1.beagle.gz -o /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/all/output -t 12 --iter 10000


awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist.all.par | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/all/indvnames.txt



cov <- as.matrix(read.table("output.cov")) # Reads estimated covariance matrix


e<-eigen(cov)

e$values/sum(e$values)


names<-read.table("/xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/par/all/indvnames.txt")

#assign the rownames of the covariance matix the 

rownames(cov)<-names$V1
#remake the plot with the colors we want

pdf(file = "pca_chr1_par_all.pdf", width = 10, height = 10, useDingbats=FALSE)

plot(e$vectors[,1:2], xlab="PC1", ylab="PC2", main="PCAngsd", col=as.factor(rownames(cov)), pch=16)

text(e$vectors[,1:2], labels=as.factor(rownames(cov)), cex= 1)

legend("topright", legend=as.factor(rownames(cov)), pch=16, col=as.factor(rownames(cov)))

dev.off()


cp /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/*/*pdf .

cp /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/pca/chr1/*/*/*pdf .


# Notes
SM1240 was included with PAR but is CRA
SM059 is missing "all" file
SM1231 was included in CRA analysis of only dups for some reason (is CRA but isn't a dup i don't think)
SM059 wasn't in the "all" analysis of PAR



# Initial Pass batch effects
# call genotypes 

#!/bin/bash

#SBATCH --job-name=initialpassgenocall
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.initialpassgenocall.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/initialpass

~/programs/angsd/angsd -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist -out /xdisk/mcnew/dannyjackson/finches/bias_testing/initialpass/genocalls -minmapQ 30 -minQ 30 -doCounts 1 -doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -SNP_pval 1e-6 -GL 1 -doMaf 2 -r NC_044571.1 -nThreads 10 

sbatch initialpass_genocall.sh 
Submitted batch job 1952897

# plink check for bias in individual missingness, frequencies, heterozygosity, and depth

plink --tped /xdisk/mcnew/dannyjackson/finches/bias_testing/initialpass/genocalls.tped --tfam /xdisk/mcnew/dannyjackson/finches/bias_testing/initialpass/genocalls.tfam --allow-extra-chr --missing --het --freq

# what do i do with frq files?

 
# age / degradation

#!/bin/bash

#SBATCH --job-name=hetero
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.hetero.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/age/

~/programs/angsd/angsd -dosaf 1 -out /xdisk/mcnew/dannyjackson/finches/bias_testing/age/chr1_noTrans.sfs -nThreads 10 -bam /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist -minQ 30 -minMapQ 30 -setMinDepthInd 4 -r NC_044571.1 -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -gl 1 -noTrans 1

~/programs/angsd/misc/realSFS chr1_noTrans.sfs.saf.idx > chr1_noTrans.est.ml

sbatch angsd_heterozygosity_age_chr1.sh 
Submitted batch job 1952671





# trying out individual bam files heterozygosity

awk 'BEGIN {FS = "/"} {print $9}' /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/bam.filelist | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/heterozygosity/indvnames.txt


#!/bin/bash

#SBATCH --job-name=heteroindv
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.heteroindv.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/heterozygosity



while read -r finch;
do 
  ~/programs/angsd/angsd -dosaf 1 -out "$finch".chr1.sfs -nThreads 10 -i /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -minQ 30 -minMapQ 30 -setMinDepthInd 4 -r NC_044571.1 -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -gl 1

  ~/programs/angsd/misc/realSFS "$finch".chr1.sfs.saf.idx > "$finch".chr1.est.ml

done < /xdisk/mcnew/dannyjackson/finches/bias_testing/dups/heterozygosity/indvnames.txt


sbatch heterozygosity.indv.sh 
Submitted batch job 1952887


a<-scan("chr1.est.ml")
a[2]/sum(a)

0.001213037

scp -r dannyjackson@filexfer.hpc.arizona.edu:/xdisk/mcnew/dannyjackson/copythis/ .



while read -r finch;
do 
  echo "all" $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt


#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.depth.%j

module load parallel
module load samtools

depthing () {
echo "depthing" "$@" >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/progress.txt 


samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$@".realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/"$@".realigned.depthstats



echo "done " $@ >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/progress.txt 
}


export -f depthing 

parallel -j 12 depthing :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  


sbatch depth.sh 
Submitted batch job 1955601

# didn't get "all" files, trying again

#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --ntasks=8
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.depth.%j

module load parallel
module load samtools

depthing () {
echo "depthing all " "$@" >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/progress.txt 


samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$@"_all.realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/"$@"_all.realigned.depthstats



echo "done all" $@ >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/progress.txt 
}


export -f depthing 

parallel -j 8 depthing :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt


sbatch depth_all.sh 
Submitted batch job 9862116

module load samtools
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/


samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1240_all.realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/SM1240_all.realigned.depthstats

samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/SM1266_all.realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/SM1266_all.realigned.depthstats













#!/bin/bash

#SBATCH --job-name=depth
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=200gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.depth.%j

while read -r finch;
do
echo "$finch" >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/depthstats.may22.txt

awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/"$finch".realigned.depthstats >> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/depthstats.may22.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplenames.txt

sbatch depthstats.sh 
Submitted batch job 9860072

# with only 50gb
9864743


# format it with tabs

sed -e :a -e '$!N;s/\nAverage = / /;ta' -e 'P;D' depthstats.may22.txt > depthstats.may22.temp.txt 

echo -e 'INDV\tAverage\tStdev' > depthstats.may22.formatted.txt

sed -e :a -e '$!N;s/\nStdev = / /;ta' -e 'P;D' depthstats.may22.temp.txt >> depthstats.may22.formatted.txt


# adding in SM059_all
awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/SM059_all.realigned.depthstats 

>> /xdisk/mcnew/dannyjackson/finches/bias_testing/depth/depthstats.may22.txt

SM059_all.realigned.depthstats