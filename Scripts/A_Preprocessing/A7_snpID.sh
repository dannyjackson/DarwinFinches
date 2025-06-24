# A6_snpID 

cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists

ls /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/indelrealignment/*bam > allsamplebams.txt


#!/bin/bash

#SBATCH --job-name=snps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.snps.%j

cd /xdisk/mcnew/finches/dannyjackson/finches/reference_lists

~/programs/angsd_elgato/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/finches/dannyjackson/finches/reference_lists/allsamplebams.txt -out /xdisk/mcnew/finches/dannyjackson/finches/reference_lists/allsnps -nThreads 10 

sbatch angsd_snps.sh 
Submitted batch job 12359997

zcat allsnps.mafs.gz | awk '{print $1, $2, $3, $4}' > sites.mafs

tail -n +2 sites.mafs > sites_headless.mafs

~/programs/angsd_elgato/angsd/angsd sites index sites_headless.mafs