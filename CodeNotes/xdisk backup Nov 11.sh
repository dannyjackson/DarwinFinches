xdisk backup Nov 11

# 20T lol
1.2T	cardinals
256K	copylist
76K	copythis
12T	finches
6.2T	finch_wgs
4.0K	pcangsd.log
2.1G	programs
4.0K	sabrina_backup
4.0K	sulidae


# go in descending order of size
#finches

12T	bias_testing
344G	vcfs # remove
295G	angsty
42G	sweed
29G	dxy
22G	vcf_likelihoods
8.3G	dadi
4.5G	reference_data
1.2G	summarystats

316M	cra
4.0K	fastqs
127M	for
280K	genelists
127M	par
207M	PCA
204K	raxml
200K	reference_lists
6.2M	timesweeper
567M	zchrom



# 12T	bias_testing
9.1T	batchnaive
2.1T	depth
# don't edit these:
18G	dups
3.5G	age
1.2G	initialpass
196M	platform
4.0K	samplesize

# batchnaive
2.5T	clipoverlap
# remove bams but keep bamstats and depthstats

#!/bin/bash

#SBATCH --job-name=zipping
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.zipping.%j

module load parallel
module load samtools

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap

ls *depthstats | parallel -j 15 bgzip

sbatch zipping.sh 
Submitted batch job 10767434


2.5T	sortedmarkedbamfiles
# remove bams but keep bamstats and depthstats
#!/bin/bash

#SBATCH --job-name=zipping
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.zipping.%j

module load parallel
module load samtools

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles

ls *depthstats | parallel -j 15 bgzip

sbatch zipping.sh 
Submitted batch job 10767376


# remove these directories
1.8T	samfiles
345G	adaptertrimmed_fastas
399G	bamfiles
352G	polygtrimmed_fastas
300G	finaltrim_fastas
269G	sortedbamfiles
197G	trim_eriks_trimmomatic_settings
128G	developing_pipeline
1.7G	indelmaps

rm -r samfiles
rm -r adaptertrimmed_fastas
rm -r bamfiles
rm -r polygtrimmed_fastas
rm -r finaltrim_fastas
rm -r sortedbamfiles
rm -r trim_eriks_trimmomatic_settings
rm -r developing_pipeline
rm -r indelmaps

# KEEP THESE
448G	indelrealignment
780K	duplicatemetrics

# depth

#!/bin/bash

#SBATCH --job-name=zipping
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=20:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.zipping.%j

module load parallel
module load samtools

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/depth

ls *depthstats | parallel -j 15 bgzip

sbatch zipping.sh 
Submitted batch job 2292602
Submitted batch job 2292607





# angsty cleaning up
238G	finches_initialcheck.tped
47G	analyses
5.2G	finches_initialcheck.geno.gz
4.0G	finches_initialcheck.mafs.gz
2.3G	oldfiles

16K	darwinfinches.arg
4.0K	darwinfinches.beagle.gz
4.0K	darwinfinches.mafs.gz
4.0K	depth
16K	finches_forpopgen2.arg
162M	finches_forpopgen2.mafs.gz
16K	finches_forpopgen.arg
508K	finches_forpopgen.mafs.gz
16K	finches_initialcheck.arg
4.0K	finches_initialcheck.tfam
80K	PCA
52K	test_depth



# copy my directory over to laurenpetrullo

#!/bin/bash

#SBATCH --job-name=copydanny
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.copydanny.%j


cd /xdisk/laurenpetrullo/dannyjackson/dannyjackson

rclone sync /xdisk/mcnew/dannyjackson/ . --progress -L --multi-thread-streams=15

sbatch copydanny.sh 
Submitted batch job 2295936



# copy isabella's directory over to laurenpetrullo

#!/bin/bash

#SBATCH --job-name=copydanny
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.copydanny.%j


cd /xdisk/laurenpetrullo/dannyjackson/iweiler

rclone sync /xdisk/mcnew/iweiler/ . --progress -L --multi-thread-streams=15

sbatch copyisabella.sh 
Submitted batch job 10779631



#!/bin/bash

#SBATCH --job-name=copylogan
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.copylogan.%j


cd /xdisk/laurenpetrullo/dannyjackson/ljvossler

rclone sync /xdisk/mcnew/ljvossler/ . --progress -L --multi-thread-streams=15

sbatch copylogan.sh 
Submitted batch job 10779631


#!/bin/bash

#SBATCH --job-name=copylogan
#SBATCH --ntasks=15
#SBATCH --nodes=1             
#SBATCH --time=120:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=100gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.copylogan.%j


cd /xdisk/laurenpetrullo/dannyjackson/ljvossler

rclone sync /xdisk/mcnew/ljvossler/ . --progress -L --multi-thread-streams=15







