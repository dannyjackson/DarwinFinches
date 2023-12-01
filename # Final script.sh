# Final script

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

# align and sort 


sbatch ~/programs/slurmscripts/alignsort.slurm

sbatch ~/programs/slurmscripts/vcfstats.slurm -i /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf


# filter VCF

bcftools view -i 'QUAL>100' /xdisk/mcnew/dannyjackson/vcfs/darwinfinches_snps_multiallelic.vcf  > darwinfinches_qualitysort.vcf

vcftools --vcf darwinfinches_qualitysort.vcf --min-meanDP 2 --remove-indels --recode --out darwinfinches_filtered


bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.recode.vcf > b10k_filtered.recode.stats.txt
bcftools stats -v /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf > b10k_filtered.geno25.maf1.stats.txt


vcftools --vcf b10k_qualitysort.vcf --min-meanDP 2 --remove-indels --recode --out b10k_filtered

sbatch ~/programs/slurmscripts/filtervcf.slurm -i /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_snps_multiallelic.vcf -n darwinfinches_nomaxDP
squeue --job 1690206

sed -i 's/\SRR2917/SRR/g' darwinfinches_filtered.recode.vcf

plink --vcf darwinfinches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.2 --maf 0.01 --recode vcf-iid --out darwinfinches_filtered_mind2


plink --vcf darwinfinches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out darwinfinches_filtered.geno25.maf1

plink --vcf darwinfinches_filtered.recode.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.02 --mind 0.1 --maf 0.01 --recode vcf-iid --out darwinfinches_filtered_cleaned

~/programs/vcf2phylip/vcf2phylip.py -i /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf


# subset by pre/post

# pre
bcftools view -s PARV1,PARV2,PL15,PL16,PL4,PL7,PL9,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/pre.vcf

# post

bcftools view -s JP4481,JP5410,JP9655,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SM1067,SM1083,SM1156,SM1157,SM1200,SM1204,SM1231,SM1237,SM1240,SM1266,SM1270,SM1271,SM1272,SM1273 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf --force-samples > /xdisk/mcnew/dannyjackson/finches/vcfs/post.vcf

# subset by species
# par
bcftools view -s PARV1,PARV2,RHC097,RHC507,SM031,SM032,SM040,SM059,SM079,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par.vcf

bcftools view -s RHC097,RHC507,SM031,SM032,SM040,SM059,SM079 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_post.vcf

bcftools view -s PARV1,PARV2,SRR329,SRR330,SRR331,SRR332,SRR333,SRR334,SRR335,SRR336,SRR337,SRR338 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/par_pre.vcf


# cra
bcftools view -s JP4481,JP5410,JP9655,PL15,PL16,PL4,PL7,PL9,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf

bcftools view -s PL15,PL16,PL4,PL7,PL9 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_pre.vcf


bcftools view -s JP4481,JP5410,JP9655,SM1067,SM1157,SM1200,SM1231,SM1240,SM1266 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/cra_post.vcf

plink --vcf /xdisk/mcnew/dannyjackson/finches/vcfs/cra.vcf --allow-extra-chr --snps-only 'just-acgt' --geno 0.25 --maf 0.1 --recode vcf-iid --indep-pairwise 50 5 0.5 --out cra_filtered.geno25.maf1



# for
bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273,SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for.vcf

bcftools view -s SRR289,SRR290,SRR291,SRR292,SRR293,SRR294,SRR295,SRR296,SRR297,SRR298 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_pre.vcf

bcftools view -s SM1083,SM1156,SM1204,SM1237,SM1270,SM1271,SM1272,SM1273 /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/for_post.vcf


# Run RAxML

cd /xdisk/mcnew/dannyjackson/finches/vcfs/

mkdir pruned 

cd pruned 

module load python

python ~/programs/sula/filter_invariants_all.py ../darwinfinches_filtered.geno25.maf1.min4.phy
mv variantsites.phy ../
mv variantsites_kept.txt ../
cd .. 
rm -r pruned

sbatch ~/programs/slurmscripts/prunevariants.slurm
squeue -j 1690221

cd /xdisk/mcnew/dannyjackson/finches/raxml

echo '10000' > p1.txt
echo '[asc~p1.txt], ASC_DNA, p1 = 1-1918469' > partitionfile.txt

# run this before running sbatch to check on p1 max

~/programs/standard-RAxML/raxmlHPC -m ASC_GTRCAT --asc-corr felsenstein -f d -d -k -n darwinfinches -q /xdisk/mcnew/dannyjackson/finches/raxml/partitionfile.txt -s /xdisk/mcnew/dannyjackson/finches/vcfs/variantsites.phy -T 12 -p 12345 -N 10 ­-b 12345 -V


sbatch ~/programs/slurmscripts/raxml.slurm



# PCA
# all
# PCA All
echo -e "CRA" > pops.txt
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

~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf  -o /xdisk/mcnew/dannyjackson/finches/PCA/all/ -p /xdisk/mcnew/dannyjackson/finches/PCA/all/pops.txt -n all -s y

cd /xdisk/mcnew/dannyjackson/finches/PCA/all/
sbatch ~/programs/slurmscripts/PCA_all.slurm

squeue -j 1690222




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


~/programs/genomics/PCA_r.sh -v /xdisk/mcnew/dannyjackson/finches/vcfs/cra_filtered.geno25.maf1.vcf -o /xdisk/mcnew/dannyjackson/finches/PCA/cra_filtered/ -p /xdisk/mcnew/dannyjackson/finches/PCA/cra/pops.txt -n cra_filtered -s y


cd /xdisk/mcnew/dannyjackson/finches/PCA/cra/
sbatch ~/programs/slurmscripts/PCA_cra.slurm
squeue --job 1682677

post
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

JP4481	JP5410	JP9655	PL15	PL16	PL4	PL7	PL9	SM1067	SM1157	SM1200	SM1231	SM1240	SM1266

JP4481  CRA post
JP5410  CRA post
JP9655  CRA post
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
squeue -j 8644180


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

PARV1	PARV2	RHC097	RHC507	SM031	SM032	SM040	SM059	SM079	SRR329	SRR330	SRR331	SRR332	SRR333	SRR334	SRR335	SRR336	SRR337	SRR338

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

squeue --job 1690225

SM1083	SM1156	SM1204	SM1237	SM1270	SM1271	SM1272	SM1273	SRR289	SRR290	SRR291	SRR292	SRR293	SRR294	SRR295	SRR296	SRR297	SRR298

SM1083  FOR post
SM1156  FOR post
SM1204  FOR post
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








# FST

grep 'NC_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf | awk '{print $1}' | sort -u | wc -l

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

grep 'NW_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf | awk '{print $1}' | sort -u | wc -l

grep -e '##' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf > /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.chroms.vcf

grep -e 'CHROM' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf >> /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.chroms.vcf


grep -e 'NW_' /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.vcf >> /xdisk/mcnew/dannyjackson/finches/vcfs/darwinfinches_filtered.geno25.maf1.chroms.vcf

# PAR

cd /xdisk/mcnew/dannyjackson/finches/reference_lists
echo -e "PARV1" > par_pre_pops.txt 
echo -e "PARV2" >> par_pre_pops.txt 
echo -e "RHC097" >> par_post_pops.txt 
echo -e "RHC507" >> par_post_pops.txt 
echo -e "SM031" >> par_post_pops.txt 
echo -e "SM032" >> par_post_pops.txt 
echo -e "SM040" >> par_post_pops.txt 
echo -e "SM059" >> par_post_pops.txt 
echo -e "SM079" >> par_post_pops.txt 
echo -e "SRR329" >> par_pre_pops.txt 
echo -e "SRR330" >> par_pre_pops.txt 
echo -e "SRR331" >> par_pre_pops.txt 
echo -e "SRR332" >> par_pre_pops.txt 
echo -e "SRR333" >> par_pre_pops.txt 
echo -e "SRR334" >> par_pre_pops.txt 
echo -e "SRR335" >> par_pre_pops.txt 
echo -e "SRR336" >> par_pre_pops.txt 
echo -e "SRR337" >> par_pre_pops.txt 
echo -e "SRR338" >> par_pre_pops.txt 



# CRA
echo -e "JP4481" > cra_post_pops.txt 
echo -e "JP5410" >> cra_post_pops.txt 
echo -e "JP9655" >> cra_post_pops.txt 
echo -e "PL15" > cra_pre_pops.txt 
echo -e "PL16" >> cra_pre_pops.txt 
echo -e "PL4" >> cra_pre_pops.txt 
echo -e "PL7" >> cra_pre_pops.txt 
echo -e "PL9" >> cra_pre_pops.txt 
echo -e "SM1067" >> cra_post_pops.txt 
echo -e "SM1157" >> cra_post_pops.txt 
echo -e "SM1200" >> cra_post_pops.txt 
echo -e "SM1231" >> cra_post_pops.txt 
echo -e "SM1240" >> cra_post_pops.txt 
echo -e "SM1266" >> cra_post_pops.txt 


# FOR
echo -e "SM1083" > for_post_pops.txt 
echo -e "SM1156" >> for_post_pops.txt 
echo -e "SM1204" >> for_post_pops.txt 
echo -e "SM1237" >> for_post_pops.txt 
echo -e "SM1270" >> for_post_pops.txt 
echo -e "SM1271" >> for_post_pops.txt 
echo -e "SM1272" >> for_post_pops.txt 
echo -e "SM1273" >> for_post_pops.txt 
echo -e "SRR289" > for_pre_pops.txt 
echo -e "SRR290" >> for_pre_pops.txt 
echo -e "SRR291" >> for_pre_pops.txt 
echo -e "SRR292" >> for_pre_pops.txt 
echo -e "SRR293" >> for_pre_pops.txt 
echo -e "SRR294" >> for_pre_pops.txt 
echo -e "SRR295" >> for_pre_pops.txt 
echo -e "SRR296" >> for_pre_pops.txt 
echo -e "SRR297" >> for_pre_pops.txt 
echo -e "SRR298" >> for_pre_pops.txt 






# Selection Scans

sbatch ~/programs/slurmscripts/selection_scans.slurm
squeue --job 8666598
