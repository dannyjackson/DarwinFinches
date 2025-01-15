# RAiSD

### Giving up on SweeD and trying RAiSD
/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -h



#!/bin/bash

#SBATCH --job-name=RAiSD_for_post
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_for_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n for_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post.vcf -M 1 -y 2 -f -p -O -R -P -D -a 1500 
# -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

# seeing if removing the -C flag helps
sbatch raised_for_post.sh 
Submitted batch job 10116808


# -M 1 imputes missing data
# -y tells it the ploidy
# -f overwrites existing files with same ID 
# -p makes a file RAiSD_Samples.STRING with sample names used in run (sanity check)
# -O shows progress in display
# -R shows additionaal info in report files
# -P makes plots
# -D generates site report
# -C path to ref fasta


/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n for_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt4.vcf  -y 2 -f -p -O -R -P -D -a 1500 

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n for_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt2.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa



# RAISD 
# FOR POST 
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


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/


# opt 3 (same as opt 2 but without regions file)
~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

#!/bin/bash

#SBATCH --job-name=RAiSD_for_post
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_for_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n for_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for_post_opt3.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa


# FOR PRE
#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/for_pre_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf

sbatch for_pre_vcf.sh 
Submitted batch job 2046340

#!/bin/bash

#SBATCH --job-name=RAiSD_for_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_for_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n for_pre -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/for/pre/for_pre.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

sbatch raisd_for_pre.sh 
Submitted batch job 2046665

# CRA PRE
#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_pre_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf

sbatch cra_pre_vcf.sh 
Submitted batch job 2046341


#!/bin/bash

#SBATCH --job-name=RAiSD_cra_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n cra_pre -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/pre/cra_pre.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa


# CRA POST
#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/cra_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf

sbatch cra_post_vcf.sh 
Submitted batch job 2046365


#!/bin/bash

#SBATCH --job-name=RAiSD_cra_post
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_cra_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n cra_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/cra/post/cra_post.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa


# PAR
# par PRE
#!/bin/bash

#SBATCH --job-name=createvcf
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.createvcf.%j


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/par_pre_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

module load bcftools 

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf

sbatch par_pre_vcf.sh 
Submitted batch job 2046342


#!/bin/bash

#SBATCH --job-name=RAiSD_par_pre
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_par_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n par_pre -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/pre/par_pre.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

sbatch raisd_par_pre.sh 
Submitted batch job 2046662


# par POST
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


cd /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post

module load bcftools 

~/programs/angsd/angsd -b /xdisk/mcnew/dannyjackson/finches/reference_lists/par_post_bams.txt -doBcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -doGeno 4 

bcftools convert /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/angsdput.bcf  -o /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf -O v 

sed -i 's/\/xdisk\/mcnew\/dannyjackson\/finches\/bias_testing\/batchnaive\/indelrealignment\///g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf

sed -i 's/\.realigned\.bam//g' /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf

sbatch par_post_vcf.sh 
Submitted batch job 2046363



#!/bin/bash

#SBATCH --job-name=RAiSD_par_post
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=300gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.RAiSD_par_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post
module load R

/home/u15/dannyjackson/programs/RAiSD/raisd-master/RAiSD -n par_post -I /xdisk/mcnew/dannyjackson/finches/sweed/vcfs/par/post/par_post.vcf -f -O -R -P -a 1500 -C /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa

sbatch raisd_par_post.sh 
Submitted batch job 2046663


# RAiSD Output
midpos start end VAR SFS LD U

while IFS=$'\t' read -r col1 col2 col3 col4 col5 col6
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
      echo "$scaf $col1 $col2 $col3 $col4 $col5 $col6" >> RAiSD_Report.par_post_manhattan
    fi
done < SweeD_Report.noca_rural




ls * | awk 'BEGIN {FS = "."} {print $3"."$4}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt


# for pre
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.for_pre.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.for_pre."$chrom" >> RAiSD_Report.for_pre.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt


# for post
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.for_post.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.for_post."$chrom" >> RAiSD_Report.for_post.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt


# cra pre
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.cra_pre.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.cra_pre."$chrom" >> RAiSD_Report.cra_pre.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt


# cra post
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.cra_post.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.cra_post."$chrom" >> RAiSD_Report.cra_post.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt

# par pre
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.par_pre.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.par_pre."$chrom" >> RAiSD_Report.par_pre.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt

# par post
cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes

echo -e "chrom\tmidpos\tstart\tend\tVAR\tSFS\tLD\tU" > RAiSD_Report.par_post.chromosomes

while read -r chrom;
do 
  awk -v chrom="$chrom" '{print chrom, $0}' RAiSD_Report.par_post."$chrom" >> RAiSD_Report.par_post.chromosomes

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/chromlist.txt

# cra pre

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes

cp RAiSD_Report.cra_pre.chromosomes RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.cra_pre.chromosomes.craplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.cra_pre.chromosomes.craplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.cra_pre.chromosomes', sep ='\t')
max(raisd$U)
# 102.8
min(raisd$U)
# 5.722e-23

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 2073013 snps, so top 0.1% would be 2073

outlier_u_disorder <- ordered_u[1:2073,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 14.01
max(outlier_u_disorder$U)
# 102.8
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.cra_pre.chromosomes.craplot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.pre.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 14.01) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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

# cra post

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes

cp RAiSD_Report.cra_post.chromosomes RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.cra_post.chromosomes.craplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.cra_post.chromosomes.craplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.cra_post.chromosomes', sep ='\t')
max(raisd$U)
# 38.41
min(raisd$U)
# 3.371e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 2141607 snps, so top 0.1% would be 2142

outlier_u_disorder <- ordered_u[1:2142,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 8.918
max(outlier_u_disorder$U)
# 38.41
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.cra_post.chromosomes.craplot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#4EAFAF", "#082B64")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions cra x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("cra.post.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 8.918) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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


# for pre

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes

cp RAiSD_Report.for_pre.chromosomes RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.for_pre.chromosomes.forplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.for_pre.chromosomes.forplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.for_pre.chromosomes', sep ='\t')
max(raisd$U)
# 68.45
min(raisd$U)
# 2.06e-22

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 5010747 snps, so top 0.1% would be 5011

outlier_u_disorder <- ordered_u[1:5011,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 8.595
max(outlier_u_disorder$U)
# 68.45
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.for_pre.chromosomes.forplot', sep ='\t')
lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.pre.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 8.595) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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


# for post

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes

cp RAiSD_Report.for_post.chromosomes RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.for_post.chromosomes.forplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.for_post.chromosomes.forplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.for_post.chromosomes', sep ='\t')
max(raisd$U)
# 45.53
min(raisd$U)
# 5.415e-14

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 5174288 snps, so top 0.1% would be 5174

outlier_u_disorder <- ordered_u[1:5174,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 8.541
max(outlier_u_disorder$U)
# 45.53
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.for_post.chromosomes.forplot', sep ='\t')

lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#FF817E", "#75002B")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("for.post.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 8.541) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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
# par pre

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes

cp RAiSD_Report.par_pre.chromosomes RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.par_pre.chromosomes.parplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.par_pre.chromosomes.parplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.par_pre.chromosomes', sep ='\t')
max(raisd$U)
# 69.78
min(raisd$U)
# 2.06e-22

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 4639528 snps, so top 0.1% would be 5011

outlier_u_disorder <- ordered_u[1:5011,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 17.61
max(outlier_u_disorder$U)
# 69.78
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.par_pre.chromosomes.parplot', sep ='\t')
lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions par x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.pre.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 17.61) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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

# par post

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes

cp RAiSD_Report.par_post.chromosomes RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044571.1 /1\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044572.1 /2\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044573.1 /3\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044574.1 /4\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044575.1 /5\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044576.1 /6\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044577.1 /7\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044578.1 /8\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044579.1 /9\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044580.1 /10\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044581.1 /11\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044582.1 /12\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044583.1 /13\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044584.1 /14\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044585.1 /15\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044586.1 /1.1\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044587.1 /17\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044588.1 /18\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044589.1 /19\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044590.1 /20\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044591.1 /21\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044592.1 /22\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044593.1 /23\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044594.1 /24\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044595.1 /25\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044596.1 /26\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044597.1 /27\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044598.1 /28\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044599.1 /29\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044600.1 /4.1\t/g' RAiSD_Report.par_post.chromosomes.parplot
sed -i 's/NC_044601.1 /Z\t/g' RAiSD_Report.par_post.chromosomes.parplot




library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
raisd <- read.csv('RAiSD_Report.par_post.chromosomes', sep ='\t')
max(raisd$U)
# 198.2
min(raisd$U)
# 1.402e-22

ordered_u <- raisd %>% 
 # desc orders from largest to smallest
 arrange(desc(U)) 

nrow(raisd)
# 4151842 snps, so top 0.1% would be 4152

outlier_u_disorder <- ordered_u[1:5011,]

outlier_u <- outlier_u_disorder %>% arrange(chrom, midpos)

min(outlier_u_disorder$U)
# 19.53
max(outlier_u_disorder$U)
# 198.2
write.csv(outlier_u, "outlieru.tsv")



# draw it with cutoff line 

# Libraries ====
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

df <- read.csv('RAiSD_Report.par_post.chromosomes.parplot', sep ='\t')
lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df$x <- as.numeric(df$chrom)
lvls <- stringr::str_sort(unique(df$chrom), numeric = TRUE)
df$chrom <- factor(df$chrom, levels = lvls)

# df <- read.csv('slidingwindow_dxy.txt', sep ='\t')
blues <- c("#A6C965", "#203000")
df.tmp <- df %>% 
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len=max(midpos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(df, ., by=c("chrom"="chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, midpos) %>%
  mutate( BPcum=midpos+tot) 
  
# get chromosome center positions par x-axis
axisdf <- df.tmp %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png("par.post.raisd.sigline.png", width=2000, height=500)

ggplot(df.tmp, aes(x=BPcum, y=(U))) +
  # Show all points
  geom_point(aes(color=as.factor(chrom)), alpha=0.8, size=.5) +
  scale_color_manual(values = rep(blues, 22 )) +

  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks= axisdf$center, guide = guide_axis(n.dodge = 2) ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,max(df$U))) + # expand=c(0,0)removes space between plot area and x axis 
  
  # add plot and axis titles
  ggtitle(NULL) +
  labs(x = "Chromosome", y = "U") +
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = 19.53) +
  # geom_hline(yintercept = -log10(sugg), linetype="dashed") +
  
  # Add highlighted points
  #geom_point(data=subset(df.tmp, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(midpos), alpha=0.7), size=5, force=1.3) +
  
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

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/for.pre.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes/for.post.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/cra.pre.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes/cra.post.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/par.pre.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes/par.post.raisd.sigline.png /xdisk/mcnew/dannyjackson/copythis


#!/bin/bash

#SBATCH --job-name=gene_raisd_cra_pre
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_cra_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.pre.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  cra.pre.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  cra.pre.relevantgenenames_15kb.txt



#!/bin/bash

#SBATCH --job-name=gene_raisd_cra_post
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_cra_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> cra.post.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  cra.post.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  cra.post.relevantgenenames_15kb.txt


#### 
#!/bin/bash

#SBATCH --job-name=gene_raisd_for_pre
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_for_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.pre.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  for.pre.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  for.pre.relevantgenenames_15kb.txt



#!/bin/bash

#SBATCH --job-name=gene_raisd_for_post
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_for_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> for.post.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  for.post.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  for.post.relevantgenenames_15kb.txt


# par pre
#!/bin/bash

#SBATCH --job-name=gene_raisd_par_pre
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_par_pre.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.pre.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  par.pre.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  par.pre.relevantgenenames_15kb.txt

# par post
#!/bin/bash

#SBATCH --job-name=gene_raisd_par_post
#SBATCH --ntasks=1
#SBATCH --nodes=1             
#SBATCH --time=5:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=50gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.gene_raisd_par_post.%j

cd /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes/

cat outlieru.tsv  | sed 's/\"//g' | sed 's/\,/\t/g' | cut -f2- > outlier_raisd.tsv

sed -i 's/ /\t/g' outlier_raisd.tsv

while read -r line;
do


chr=`awk 'BEGIN {FS = "\t"} {print $1}' <<<"${line}"`
midpos=`awk 'BEGIN {FS = "\t"} {print $2}' <<<"${line}"`
minpos=$((midpos - 15000))
maxpos=$((midpos + 15000))

grep "$chr" /xdisk/mcnew/dannyjackson/finches/reference_data/ncbi_dataset/data/GCF_901933205.1/genomic.gff | grep 'ID=gene' | awk '{OFS = "\t"} ($4 < '$maxpos' && $4 > '$minpos' || $5 < '$maxpos' && $5 > '$minpos' || $4 < '$minpos' && $5 > '$minpos' || $4 < '$maxpos' && $5 > '$maxpos')' >> par.post.relevantgenes_15kb.txt


done < <(tail -n +2 outlier_raisd.tsv)




awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}'  par.post.relevantgenes_15kb.txt | sed 's/ID\=gene\-//g' | sort -u >  par.post.relevantgenenames_15kb.txt




sbatch genes_cra_post.sh 
sbatch genes_cra_pre.sh  
sbatch gene_for_post.sh  
sbatch gene_for_pre,sh  
sbatch gene_par_post.sh  
sbatch gene_par_pre.sh   

sbatch gene_for_post.sh  
Submitted batch job 10137054
sbatch gene_for_pre.sh  
Submitted batch job 10137068
batch gene_par_post.sh  
Submitted batch job 10137055
sbatch gene_par_pre.sh   
Submitted batch job 10137057



cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/cra.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes/cra.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/for.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes/for.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/par.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes/par.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/copythis

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/cra.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes/cra.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/for.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes/for.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/par.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists
cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes/par.post.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists


wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/cra.pre.relevantgenenames_15kb.txt 
wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/post/reports/chromosomes/cra.post.relevantgenenames_15kb.txt 
wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/for.pre.relevantgenenames_15kb.txt
wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/post/reports/chromosomes/for.post.relevantgenenames_15kb.txt 
wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/par.pre.relevantgenenames_15kb.txt 
wc -l /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/post/reports/chromosomes/par.post.relevantgenenames_15kb.txt 


cp -r /xdisk/mcnew/dannyjackson/finches/genelists /xdisk/mcnew/dannyjackson/copylist


cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/cra/pre/reports/chromosomes/cra.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/cra.pre.relevantgenenames_RAiSD.txt

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/for/pre/reports/chromosomes/for.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/for.pre.relevantgenenames_RAiSD.txt

cp /xdisk/mcnew/dannyjackson/finches/sweed/RAiSD/par/pre/reports/chromosomes/par.pre.relevantgenenames_15kb.txt /xdisk/mcnew/dannyjackson/finches/genelists/par.pre.relevantgenenames_RAiSD.txt
