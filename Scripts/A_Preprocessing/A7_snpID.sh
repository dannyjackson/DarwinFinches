# SNP ID
This uses ANGSD to identify SNPs from the bam files developed in the previous analyses.
---

First, create list of reference bams
```bash
cd /xdisk/mcnew/dannyjackson/cardinals_dfinch/referencelists

ls /xdisk/mcnew/dannyjackson/cardinals_dfinch/datafiles/indelrealignment/*bam > allsamplebams.txt

```
Create sbatch script
```bash
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

cd /xdisk/mcnew/finches/dannyjackson/finches/referencelists

~/programs/angsd_elgato/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/finches/dannyjackson/finches/reference_lists/allsamplebams.txt -out /xdisk/mcnew/finches/dannyjackson/finches/reference_lists/allsnps -nThreads 10 

```
Submit
```bash
sbatch angsd_snps.sh 
```
Filter to a sites file that only retains relevant columns, no header,and is indexed.
```bash
zcat allsnps.mafs.gz | awk '{print $1, $2, $3, $4}' > sites.mafs

tail -n +2 sites.mafs > sites_headless.mafs

~/programs/angsd_elgato/angsd/angsd sites index sites_headless.mafs
```

# Compute putative polymorphic sites in each species

#!/bin/bash

#SBATCH --job-name=crasnps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.crasnps.%j

cd /xdisk/mcnew/finches/dannyjackson/finches/referencelists

~/programs/angsd_elgato/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 7 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/cra_all_bams.txt -out /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crasnps -nThreads 10 -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs

#!/bin/bash

#SBATCH --job-name=forsnps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.forsnps.%j

cd /xdisk/mcnew/finches/dannyjackson/finches/referencelists

~/programs/angsd_elgato/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 9 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_all_bams.txt -out /xdisk/mcnew/finches/dannyjackson/finches/referencelists/forsnps -nThreads 10 -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs

#!/bin/bash

#SBATCH --job-name=parsnps
#SBATCH --ntasks=10
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.parsnps.%j

cd /xdisk/mcnew/finches/dannyjackson/finches/referencelists

~/programs/angsd_elgato/angsd/angsd -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 32 -SNP_pval 1e-6 -setMinDepthInd 4 -minInd 9 -minQ 30 -minMaf 0.05 -minMapQ 30 -bam /xdisk/mcnew/finches/dannyjackson/finches/referencelists/par_all_bams.txt -out /xdisk/mcnew/finches/dannyjackson/finches/referencelists/parsnps -nThreads 10 -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs
