# Sweed

#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j

# ~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 2 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -out /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa -nThreads 8 -doBcf 1 -doPost 1

cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs

# opt 1
# ~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs

# base file
# ~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 3 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs

# opt 3 (same as opt 2 but without regions file)
~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

sbatch for_post_opt3.sh 
Submitted batch job 2046290

mv angsdput.bcf for_post_opt3.bcf

mv angsdput.mafs.gz for_post_opt3.mafs.gz

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.vcf

# opt 2
~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 -r NC_044571.1

# opt 4
~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 -rf regionsfile.txt

mv angsdput.bcf for_post_opt4.bcf


bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt4.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt4.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt4.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt4.vcf

NC_044571.1
NC_044572.1

# opt 3
sbatch for_post_createvcf.sh 
Submitted batch job 2045695

sbatch createvcf.sh 
Submitted batch job 2045385


sbatch createvcf_opt1.sh 
Submitted batch job 10112221

mv angsdput.bcf for_post_opt2.bcf

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf

~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -grid 100000 -length 100000

grep -v /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf | wc -l

# option 2
mv angsdput.bcf for_post_opt2.bcf

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf

~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -grid 15

# option 1
mv angsdput.bcf for_post_opt1.bcf

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt1.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt1.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt1.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt1.vcf

~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt1.vcf -grid 15

~/programs/SweeD/SweeD -name testing -input darwinfinches_vcftools.DP4.vcf.recode.vcf -grid 15


darwinfinches_quality_DP6.vcf
# first output 
mv angsdput.bcf for_post.bcf

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf

~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf -grid 100 

~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf -grid 15

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






#!/bin/bash

#SBATCH --job-name=sweed_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/pre


~/programs/SweeD/SweeD -name cra_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000

sbatch sweed_cra_pre.sh 
Submitted batch job 2046735



# SWEED
~/programs/SweeD/SweeD -name for_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -grid 100000 -length 100000



# CRA
#!/bin/bash

#SBATCH --job-name=sweed_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/pre


~/programs/SweeD/SweeD -name cra_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000

zcat /scratch/dnjacks4/cardinalis/to_b10k/sweed/urban_noca.filtered.geno25.maf1.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq > scaff_names_noca_urban.txt


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


#!/bin/bash

#SBATCH --job-name=sweed_cra_post
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/post


~/programs/SweeD/SweeD -name cra_post -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000

sbatch sweed_cra_post.sh 
Submitted batch job 10135259

# FOR

#!/bin/bash

#SBATCH --job-name=sweed_for_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_for_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/for/pre


~/programs/SweeD/SweeD -name for_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf -grid 100 -length 100000

sbatch for_pre.sh 
Submitted batch job 2047360

#!/bin/bash

#SBATCH --job-name=sweed_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/pre


~/programs/SweeD/SweeD -name cra_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000

# PAR

#!/bin/bash

#SBATCH --job-name=sweed_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/pre


~/programs/SweeD/SweeD -name cra_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000

#!/bin/bash

#SBATCH --job-name=sweed_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/cra/pre


~/programs/SweeD/SweeD -name cra_pre -input /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -grid 100000 -length 100000













# Redoing sweed but by chromosome length to set grid size that pins 15kb apart

# CRA PRE CHROMS
#!/bin/bash

#SBATCH --job-name=sweed_cra_pre_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf > cra.$chr.vcf


~/programs/SweeD/SweeD -name cra_pre_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/pre/cra.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch cra_pre_sweed_chroms.sh 
Submitted batch job 2047397


# CRA POST CHROMS
#!/bin/bash

#SBATCH --job-name=sweed_cra_post_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_post_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf > cra.$chr.vcf


~/programs/SweeD/SweeD -name cra_pre_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/post/cra.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch cra_post_sweed_chroms.sh 
Submitted batch job 2047399

# FOR PRE CHROMS
#!/bin/bash

#SBATCH --job-name=sweed_for_pre_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_for_pre_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf > for.$chr.vcf


~/programs/SweeD/SweeD -name for_pre_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/pre/for.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch for_pre_sweed_chroms.sh 
Submitted batch job 2047400


# FOR POST CHROMS
cp /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.vcf /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/post/for_post.vcf

#!/bin/bash

#SBATCH --job-name=sweed_for_post_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_for_post_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/post/for_post.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/post/for_post.vcf > for.$chr.vcf


~/programs/SweeD/SweeD -name for_pre_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/post/for.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch for_post_sweed_chroms.sh 
Submitted batch job 2047405

# PAR PRE CHROMS
#!/bin/bash

#SBATCH --job-name=sweed_par_pre_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_par_pre_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf > par.$chr.vcf


~/programs/SweeD/SweeD -name par_pre_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/pre/par.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch par_pre_sweed_chroms.sh 
Submitted batch job 2047402


# PAR POST CHROMS
#!/bin/bash

#SBATCH --job-name=sweed_par_post_chroms
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_par_post_chroms.%j

grep 'length' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf | grep 'NC_' > chroms.txt

sed -i 's/##contig=<ID=//g' chroms.txt
sed -i 's/,length=/\t/g' chroms.txt
sed -i 's/>/\t/g' chroms.txt

while read -r line;
do

chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
chrlength=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
gridvar=$(($chrlength / 15000))


grep "^#\|^$chr" /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf > par.$chr.vcf


~/programs/SweeD/SweeD -name par_post_$chr -input /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/post/par.$chr.vcf -grid "$gridvar" 

done < chroms.txt

sbatch par_post_sweed_chroms.sh 
Submitted batch job 2047403




# make single file in each directory that compiles all raisd output

# cra pre
#!/bin/bash

#SBATCH --job-name=sweed_cra_pre_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=30
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_pre_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/pre

grep 'Position' SweeD_Report.cra_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.cra_pre

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.cra_pre
    fi
done < SweeD_Report.cra_pre_"$chrom" 
done < chroms.txt

sbatch concat.sh 
Submitted batch job 2048751


# cra post

#!/bin/bash

#SBATCH --job-name=sweed_cra_post_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_cra_post_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/post

grep 'Position' SweeD_Report.cra_post_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.cra_post

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.cra_post
    fi
done < SweeD_Report.cra_pre_"$chrom" 
done < chroms.txt

sbatch concat.sh 
Submitted batch job 2048762

# for pre

#!/bin/bash

#SBATCH --job-name=sweed_for_pre_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_for_pre_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/pre

grep 'Position' SweeD_Report.for_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.for_pre

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.for_pre
    fi
done < SweeD_Report.for_pre_"$chrom" 
done < chroms.txt

sbatch concat.sh 
Submitted batch job 10137098


# for post

#!/bin/bash

#SBATCH --job-name=sweed_for_post_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_for_post_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/post

grep 'Position' SweeD_Report.for_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.for_post

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.for_post
    fi
done < SweeD_Report.for_pre_"$chrom" 
done < chroms.txt

sbatch concat.sh 
Submitted batch job 10137100

# par pre
#!/bin/bash

#SBATCH --job-name=sweed_par_pre_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_par_pre_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/pre

grep 'Position' SweeD_Report.par_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.par_pre

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.par_pre
    fi
done < SweeD_Report.par_pre_"$chrom" 
done < chroms.txt

sbatch concat.sh 

# par post

#!/bin/bash

#SBATCH --job-name=sweed_par_post_concat
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.sweed_par_post_concat.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/post

grep 'Position' SweeD_Report.par_post_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.par_post

while read -r line1;
do 

chrom=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line1}"`

while read -r line;
do 
col1=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`

    if [[ $col1 = 'Position' ]]; then  
        continue
    elif
      [[ $col1 = //* ]]; then
      continue
    elif
      [ -z "${col1}" ]; then
      continue
    else
      echo "$chrom $line " >> SweeD_Report.par_post
    fi
done < SweeD_Report.par_post_"$chrom" 
done < chroms.txt

sbatch concat.sh 
Submitted batch job 10137114












# find outlier regions
# CRA PRE
cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/pre

# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.cra_pre
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.cra_pre

library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.cra_pre', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 1200.883
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 158.8352
max(outlier_sweed_disorder$Likelihood)
# 1200.883
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.cra_pre SweeD_Report.cra_pre.formanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.cra_pre.formanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.cra_pre.formanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.cra_pre.formanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.pre.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 158.8352, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()


cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.pre.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp cra.pre.sweed.* /xdisk/mcnew/dannyjackson/copythis/







# CRA POST 
cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/post


# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.cra_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.cra_post

grep 'Position' SweeD_Report.cra_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.cra_post_temp

cat SweeD_Report.cra_post_temp > tmp.txt
cat SweeD_Report.cra_post >> tmp.txt


mv tmp.txt SweeD_Report.cra_post

sed -i 's/.1 /.1\t/g' SweeD_Report.cra_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.cra_post


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.cra_post', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 131.294
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 30.70069
max(outlier_sweed_disorder$Likelihood)
# 131.294
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.cra_post SweeD_Report.cra_post.formanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.cra_post.formanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.cra_post.formanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.cra_post.formanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.pre.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 30.70069, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.post.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp cra.post.sweed.* /xdisk/mcnew/dannyjackson/copythis/






# FOR PRE
cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/pre


# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.for_pre
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.for_pre

grep 'Position' SweeD_Report.for_pre_NC_044571.1 | head -n 1 | awk '{print "Chromosome", $0}' > SweeD_Report.for_pre_temp


sed -i 's/.1 /.1\t/g' SweeD_Report.for_pre
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.for_pre


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.for_pre', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 176.1513
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 33.98217
max(outlier_sweed_disorder$Likelihood)
# 176.1513
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.for_pre SweeD_Report.for_pre.formanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.for_pre.formanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.for_pre.formanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.for_pre.formanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.pre.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 33.98217, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.pre.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp for.pre.sweed.* /xdisk/mcnew/dannyjackson/copythis/


# FOR POST

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/post


# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.for_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.for_post



sed -i 's/.1 /.1\t/g' SweeD_Report.for_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.for_post


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.for_post', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 115.4548
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 29.63933
max(outlier_sweed_disorder$Likelihood)
# 115.4548
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.for_post SweeD_Report.for_post.formanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.for_post.formanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.for_post.formanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.for_post.formanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.post.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 33.98217, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.post.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp for.post.sweed.* /xdisk/mcnew/dannyjackson/copythis/



# PAR PRE

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/pre


# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.par_pre
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.par_pre



library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.par_pre', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 313.7581
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 61.62162
max(outlier_sweed_disorder$Likelihood)
# 313.7581
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.par_pre SweeD_Report.par_pre.parmanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.par_pre.parmanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.par_pre.parmanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.par_pre.parmanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions par x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.pre.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 61.62162, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp par.pre.sweed.* /xdisk/mcnew/dannyjackson/copythis/



# PAR POST

cd /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/post


# alternative subest of just the top 0.1%
sed -i 's/.1 /.1\t/g' SweeD_Report.par_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.par_post



sed -i 's/.1 /.1\t/g' SweeD_Report.par_post
sed -i 's/Chromosome /Chromosome\t/g' SweeD_Report.par_post


library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
sweed <- read.csv('SweeD_Report.par_post', sep ='\t', row.names=NULL)
max(sweed$Likelihood, na.rm=TRUE)
# 1105.309
min(sweed$Likelihood, na.rm=TRUE)
# 0

ordered_sweed <- sweed %>% 
 # desc orders from largest to smallest
 arrange(desc(Likelihood)) 

nrow(sweed)
# 67775 snps, so top 0.1% would be 68

outlier_sweed_disorder <- ordered_sweed[1:68,]

outlier_sweed <- outlier_sweed_disorder %>% arrange(Chromosome, Position)

min(outlier_sweed_disorder$Likelihood)
# 83.9322
max(outlier_sweed_disorder$Likelihood)
# 1105.309
write.csv(outlier_sweed, "outliersweed.tsv")


# plot
cp SweeD_Report.par_post SweeD_Report.par_post.parmanhattan

sed -i 's/NC_044571.1/1\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044572.1/2\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044573.1/3\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044574.1/4\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044575.1/5\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044576.1/6\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044577.1/7\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044578.1/8\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044579.1/9\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044580.1/10\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044581.1/11\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044582.1/12\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044583.1/13\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044584.1/14\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044585.1/15\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044586.1/1.1\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044587.1/17\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044588.1/18\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044589.1/19\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044590.1/20\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044591.1/21\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044592.1/22\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044593.1/23\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044594.1/24\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044595.1/25\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044596.1/26\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044597.1/27\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044598.1/28\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044599.1/29\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044600.1/4.1\t/g' SweeD_Report.par_post.parmanhattan
sed -i 's/NC_044601.1/999\t/g' SweeD_Report.par_post.parmanhattan

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.table('SweeD_Report.par_post.parmanhattan', sep ='\t', row.names=NULL, header = TRUE)

blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot) 
  
# get chromosome center positions par x-axis
axisdf <- df.tmp %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.post.sweed.15kb.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(Likelihood))) +
  # Show all points
  geom_point(aes(color=as.factor(Chromosome)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chromosome, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1)) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "CLR") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 83.9322, linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(Position), alpha=0.7), size=5, force=1.3) +
  
  # Custom the theme:
  theme_bw(base_size = 22) +
  theme( 
    plot.title = element_text(hjust = 0.5),
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


dev.off()





cat outliersweed.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_sweed.tsv

sed -i 's/ /\t/g' outlier_sweed.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
winstart=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((winstart - 15000))
maxpos=$((winstart + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.post.sweed.auto.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_sweed.tsv)

cp par.post.sweed.* /xdisk/mcnew/dannyjackson/copythis/



cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/pre/cra.pre.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/cra/post/cra.post.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists

cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/pre/for.pre.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/for/post/for.post.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists

cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/pre/par.pre.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/SweeD/chromosomes/par/post/par.post.sweed.auto.relevantgenes_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  cra.post.sweed.auto.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  cra.post.relevantgenenames_15kb.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.post.sweed.auto.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  for.post.relevantgenenames_15kb.txt

awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' for.pre.sweed.auto.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  for.pre.relevantgenenames_15kb.txt

cp for.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis/for.pre.sweed.txt



