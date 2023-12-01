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

cd /xdisk/mcnew/dannyjackson/finches/fst/par
rm *
echo -e "PARV1" > pre.txt 
echo -e "PARV2" >> pre.txt 
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
cd /xdisk/mcnew/dannyjackson/finches/fst/cra

echo -e "JP4481" > post.txt 
echo -e "JP5410" >> post.txt 
echo -e "JP9655" >> post.txt 
echo -e "PL15" > pre.txt 
echo -e "PL16" >> pre.txt 
echo -e "PL4" >> pre.txt 
echo -e "PL7" >> pre.txt 
echo -e "PL9" >> pre.txt 
echo -e "SM1067" >> post.txt 
echo -e "SM1157" >> post.txt 
echo -e "SM1200" >> post.txt 
echo -e "SM1231" >> post.txt 
echo -e "SM1240" >> post.txt 
echo -e "SM1266" >> post.txt 


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


# Nucleotide Diversity 

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


awk '$5>= 2217320 && $4>= 2203970' test.gff

# the above works! I need to incorporate it into the other awk script now ugh


/xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff

# confusing notes


awk '{$2 > 140268 && $2 < 142418}1' /scratch/dnjacks4/cardinalis/to_b10k/b10k_filtered.geno25.maf1.vcf | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01020886/' >> Col6a1.vcf



awk '{$5 > 4400001 && $4 < 4410000}1' /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | awk '$0 !~ /\##/' | awk '$0 ~ /VYXE01006626/' >> rho.vcf




sbatch ~/programs/slurmscripts/selection_scans.slurm
