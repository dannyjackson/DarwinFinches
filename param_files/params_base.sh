module load R/4.4.0
module load htslib/1.19.1
# module load bedtools2/2.29.2
module load python/3.11/3.11.4
module load bwa/0.7.18
module load bcftools/1.19
module load vcftools/0.1.16
module load plink/1.9
module load samtools/1.19.2

# Define variables
# all
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches/    # main directory for output files
PROGDIR=~/programs  # path to directory for all installed programs
BAMDIR=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/indelrealignment/
PROJHUB=DarwinFinches
SCRIPTDIR=${PROGDIR}/${PROJHUB}
PATH=$PATH:$SCRIPTDIR # this adds the workshop script directory to our path, so that executable scripts in it can be called without using the full path
ID=finches
FILENAME_LIST="${OUTDIR}/referencelists/samplenames.txt" # list with sample codes associated with each file in dataset, one per line

# define aspects of the reference genome
CHRLEAD=NC_0 # characters at the start of a chromosome number (excluding scaffolds)
SEXCHR=NC_044601
REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna # path to reference genome
GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff # path to gff file


# define the path for the chromosome conversion file (converts chromosome ascension names to numbers)
CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

# define variables for filtering in angsd
MINMAPQ=30
MINQ=30
MININD=4

source ~/programs/DarwinFinches/base_setup.sh
