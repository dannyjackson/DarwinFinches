module load bwa/0.7.17
module load samtools/1.19.2
module load bowtie2
module load picard
module load samtools
module load parallel
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9

source ~/programs/DarwinFinches/param_files/params_base.sh
source ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/preprocessing_setup.sh

# across all preprocessing
THREADS=12

# trimming
FASTAS=/path/to/fasta/files # format must be samplename_
TRIMJAR=/path/to/trimmomatic/jarfile.jar
LEAD=<SET_VALUE> # value to trim from leading strand, often 20
TRAIN=<SET_VALUE> # value to trim from trailing strand, often 20
SLIDE=<SET_VALUE> # threshold and windlow length, often 4:20
MINREADLEN=<SET_VALUE> # minimum length for a read to be kept, often 90 for 150bp sequencing

# clipping
BAMUTILBAM=/path/to/bamutil/bin/bam/file

# bam statistics

# snp ID
ANGSD=~/programs/angsd/ # path to directory with angsd executables
SNPPVAL=1e-6 # max p-value for snp to be considered significant, often 1e-6
MINDEPTHIND=4 # minimum depth per individual required for a site to be kept
MININD=4 # minimum number of individuals required for a site to be kept
MINQ=30 # minimum quality score required for a site to be kept
MINMAF=0.05 # minimum minor allele frequency required for a site to be kept
MINMAPQ=30 # minimum mapping quality score required for a site to be kept

