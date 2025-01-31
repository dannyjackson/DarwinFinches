# dadi


# make vcfs for dadi

#!/bin/bash

#SBATCH --job-name=post_bcf
#SBATCH --ntasks=8
#SBATCH --nodes=1             
#SBATCH --time=40:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.post_bcf.%j

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/postbams.txt -out /xdisk/mcnew/dannyjackson/finches/dadi/post/all/all -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1

sbatch post_bcfs.sh 
Submitted batch job 2147197



bcftools convert all.bcf  -o all.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' all.vcf

sed -i 's/\.realigned\.bam//g' all.vcf


~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx > par_post.ml

cd ../cra
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx > cra_post.ml

cd ../for
~/programs/angsd/misc/realSFS /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx > for_post.ml


~/programs/angsd/misc/realSFS dadi /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/cra/post/cra_post.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/for/post/for_post.saf.idx /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/par/post/par_post.saf.idx -sfs /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/populations/post/cra/cra_post.ml /xdisk/mcnew/dannyjackson/finches/dadi/post/for/for_post.ml /xdisk/mcnew/dannyjackson/finches/dadi/post/par/par_post.ml -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

-ref hg19NoChr.fa -anc hg19ancNoChr.fa.gz 

# generate small test file

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/par_post_bams.txt -out par -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -r NC_044595.1

~/programs/angsd/misc/realSFS par.saf.idx > par.ml

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_bams.txt -out cra -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -r NC_044595.1

~/programs/angsd/misc/realSFS cra.saf.idx > cra.ml

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -out for -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -r NC_044595.1

~/programs/angsd/misc/realSFS for.saf.idx > for.ml


~/programs/angsd/misc/realSFS dadi par.saf.idx cra.saf.idx for.saf.idx > three.sfs

cat three.sfs

~/programs/angsd/misc/realSFS dadi par.saf.idx cra.saf.idx for.saf.idx -sfs par.idx.ml -sfs cra.idx.ml -sfs for.idx.ml -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > three.dadi.sfs

echo '10 8 9 unfolded' > three.head.sfs
cat three.sfs >> three.head.sfs

~/programs/angsd/misc/realSFS par.saf.idx cra.saf.idx for.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > three.sfs

~/programs/angsd/misc/realSFS par.saf.idx cra.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > two.sfs

echo '10 8 unfolded' > two.head.sfs
cat two.sfs >> two.head.sfs

# generate total genome file

# first, make autosome regions file 
echo -e 'NC_044571.1\nNC_044572.1\nNC_044573.1\nNC_044574.1\nNC_044575.1\nNC_044576.1\nNC_044577.1\nNC_044578.1\nNC_044579.1\nNC_044580.1\nNC_044581.1\nNC_044582.1\nNC_044583.1\nNC_044584.1\nNC_044585.1\nNC_044586.1\nNC_044587.1\nNC_044588.1\nNC_044589.1\nNC_044590.1\nNC_044591.1\nNC_044592.1\nNC_044593.1\nNC_044594.1\nNC_044595.1\nNC_044596.1\nNC_044597.1\nNC_044598.1\nNC_044599.1\nNC_044600.1' > /xdisk/mcnew/dannyjackson/finches/reference_lists/autosomes.txt


#!/bin/bash

#SBATCH --job-name=post_sfs
#SBATCH --ntasks=8
#SBATCH --nodes=1             
#SBATCH --time=40:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=10gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.post_sfs.%j

cd /xdisk/mcnew/dannyjackson/finches/dadi/post/all

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/par_post_bams.txt -out par -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -rf /xdisk/mcnew/dannyjackson/finches/reference_lists/autosomes.txt

~/programs/angsd/misc/realSFS par.saf.idx > par.ml

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_bams.txt -out cra -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -rf /xdisk/mcnew/dannyjackson/finches/reference_lists/autosomes.txt

~/programs/angsd/misc/realSFS cra.saf.idx > cra.ml

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -out for -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1 -rf /xdisk/mcnew/dannyjackson/finches/reference_lists/autosomes.txt

~/programs/angsd/misc/realSFS for.saf.idx > for.ml

~/programs/angsd/misc/realSFS par.saf.idx cra.saf.idx for.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > three.sfs

echo '10 8 9 unfolded' > three.head.sfs
cat three.sfs >> three.head.sfs



# run these as tests

# species = sample size | haploid | diploid
# par = 7 | 8 | 15
# cra = 9 | 10 | 19
# for = 8 | 9 | 17

# par cra
~/programs/angsd/misc/realSFS par.saf.idx cra.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > two.sfs

echo '10 8 unfolded' > two.head.sfs
cat two.sfs >> two.head.sfs

mv two.sfs two.par_cra.sfs
mv two.head.sfs two.par_cra.head.sfs

# par for
~/programs/angsd/misc/realSFS par.saf.idx for.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > two.par_for.sfs

echo '10 9 unfolded' > two.par_for.head.sfs
cat two.par_for.sfs >> two.par_for.head.sfs


# cra for
~/programs/angsd/misc/realSFS cra.saf.idx for.saf.idx -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa > two.cra_for.sfs

echo '10 9 unfolded' > two.cra_for.head.sfs
cat two.cra_for.sfs >> two.cra_for.head.sfs


sed -i 's/folded/unfolded/g' two.cra_for.head.sfs



two.cra_for.head.sfs



# my edits
cp two.par_cra.head.sfs two.par_cra.head.chr.sfs
cp two.cra_for.head.sfs two.cra_for.head.chr.sfs
cp two.par_for.head.sfs two.par_for.head.chr.sfs


nano two.par_cra.head.sfs
nano two.par_cra.head.chr.sfs
nano two.cra_for.head.sfs 
nano two.cra_for.head.chr.sfs
nano two.par_for.head.sfs 
nano two.par_for.head.chr.sfs


import dadi
import matplotlib.pyplot as plt

fs = dadi.Spectrum.from_file('two.par_cra.head.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('PAR')
plt.xlabel('CRA')
plt.savefig('par_cra_2d_spectrum.png')

plt.clf()

fs = dadi.Spectrum.from_file('two.cra_for.head.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('CRA')
plt.xlabel('FOR')
plt.savefig('cra_for_2d_spectrum.png')


plt.clf()

fs = dadi.Spectrum.from_file('two.par_for.head.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('PAR')
plt.xlabel('FOR')
plt.savefig('par_for_2d_spectrum.png')

plt.clf()

# chromosomes

fs = dadi.Spectrum.from_file('two.par_cra.head.chr.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('PAR')
plt.xlabel('CRA')
plt.savefig('par_cra_2d_chr_spectrum.png')

plt.clf()

fs = dadi.Spectrum.from_file('two.cra_for.head.chr.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('CRA')
plt.xlabel('FOR')
plt.savefig('cra_for_2d_chr_spectrum.png')


plt.clf()

fs = dadi.Spectrum.from_file('two.par_for.head.chr.sfs')
plot_spectrum = dadi.Plotting.plot_single_2d_sfs(fs, vmin=1)
plt.ylabel('PAR')
plt.xlabel('FOR')
plt.savefig('par_for_2d_chr_spectrum.png')

plt.clf()
