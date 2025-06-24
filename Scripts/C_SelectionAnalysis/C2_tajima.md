# C4_tajima.md
# zip some files up
sbatch --account=mcnew \
        --job-name=zipping \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.zipping.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=30:00:00 \
        --mem=50gb \
        /xdisk/mcnew/dannyjackson/cardinals/analyses/thetas/zipthemfilesup.sh

# run once per population to generate all files
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )

# Iterate over each combination
for sp in "${species[@]}"; do
sbatch --account=mcnew \
        --job-name=tajima_500000_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_500000_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=10:00:00 \
        --mem=50gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
        -w 500000 -s 500000
done




# Define species, environments, and window sizes
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )
window_sizes=( 50000 5000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        step=$(expr $win / 4)

        sbatch --account=mcnew \
               --job-name=tajima_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.tajima_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=1:00:00 \
               --mem=50gb \
               ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
               -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
               -w $win -s $step
    done
done


# do snps
# Define species, environments, and window sizes
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )

# Iterate over each combination
for sp in "${species[@]}"; do

sbatch --account=mcnew \
        --job-name=tajima_${win}_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_${win}_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=1:00:00 \
        --mem=50gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
        -w 1 -s 1
done



# calculate change in Tajima
# 50,000 snps

# cra
(echo -e "chromo\tposition\tTajima"; awk 'BEGIN {OFS="\t"} NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' crapre/crapre.Tajima.50000.Ztransformed.csv crapost/crapost.Tajima.50000.Ztransformed.csv) | grep -v 'NW'  | grep -v 'Z' >> cradiff.txt


# for
(echo -e "chromo\tposition\tTajima"; awk 'NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' forpre/forpre.Tajima.50000.Ztransformed.csv forpost/forpost.Tajima.50000.Ztransformed.csv) | head >> fordiff.txt

# par
(echo -e "chromo\tposition\tTajima"; awk 'NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' parpre/parpre.Tajima.50000.Ztransformed.csv parpost/parpost.Tajima.50000.Ztransformed.csv) >> pardiff.txt

# plot these changes
# cra




# Repeat outlier and plotting analyses using depth and map filtered windows

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas

species=( "cra" "for" "par")
# first, add header to input file (double check if necessary before doing so)
for sp in "${species[@]}"; do

    PREFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}pre/50000/${sp}pre.theta.thetasWindow.pestPG.depthmapfiltered"
    
    POSTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}post/50000/${sp}post.theta.thetasWindow.pestPG.depthmapfiltered"

    (echo -e "chromo\tposition\tTajima"; awk 'BEGIN {OFS="\t"} NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' ${PREFILE} ${POSTFILE}) >> ${sp}_tajimadiff.depthmapfiltered.txt

done


species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    source ~/programs/DarwinFinches/param_files/${sp}pre_params_tajima.sh

    FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}_tajimadiff.depthmapfiltered.txt"

    cp ${FILE} "${FILE}.numchrom" 

    # Replace chromosome names if conversion file is provided
    if [ -n "$CHR_FILE" ]; then
        echo "Replacing chromosome names based on conversion file..."
        while IFS=',' read -r first second; do
            sed -i "s/$second/$first/g" "${FILE}.numchrom" 
        done < "$CHR_FILE"
    fi

    # z transform windowed data
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.filteredfiles.r" \
    "${OUTDIR}" "${CUTOFF}" "${FILE}.numchrom" 

    cp ${FILE} "${FILE}.numchrom" 

    Z_OUT="${FILE}.numchrom.Ztransformed.csv"

    # sed -i 's/\"//g' ${Z_OUT}

    # Run R script for plotting
    echo "Generating Manhattan plot from ${Z_OUT}..."
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.filteredfiles.tajimadiff.r" \
        "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "Tajima" 

    echo "Script completed successfully!"

done


# test for significant difference in Tajima's D using 50kb windows
## CRA

library(dplyr)
library(ggplot2)
library(ggfortify)

### load data
cra_pre <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/50000/crapre.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(cra_pre) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

cra_post <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/50000/crapost.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(cra_post) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

### Paired:
t.test(cra_pre$Tajima, cra_post$Tajima, paired=TRUE)
> t.test(cra_pre$Tajima, cra_post$Tajima, paired=TRUE)

        Paired t-test

data:  cra_pre$Tajima and cra_post$Tajima
t = -119.44, df = 75560, p-value < 2.2e-16
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
 -0.3431203 -0.3320407
sample estimates:
mean difference 
     -0.3375805 


## FOR

library(dplyr)
library(ggplot2)
library(ggfortify)

### load data
for_pre <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/50000/forpre.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(for_pre) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

for_post <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/50000/forpost.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(for_post) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

### Paired:
t.test(for_pre$Tajima, for_post$Tajima, paired=TRUE)

> t.test(for_pre$Tajima, for_post$Tajima, paired=TRUE)

        Paired t-test

data:  for_pre$Tajima and for_post$Tajima
t = 112.69, df = 75560, p-value < 2.2e-16
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
 0.1730137 0.1791389
sample estimates:
mean difference 
      0.1760763 
#### mean for pre = 1.062328
#### mean for post = 0.8862516

## PAR

library(dplyr)
library(ggplot2)
library(ggfortify)

### load data
par_pre <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/50000/parpre.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(par_pre) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

par_post <- read.delim("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/50000/parpost.theta.thetasWindow.pestPG.autosomes", header=FALSE) %>% na.omit()
colnames(par_post) <- c("info", "Chr", "WinCenter", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites")

### Paired:
t.test(par_pre$Tajima, par_post$Tajima, paired=TRUE)

        Paired t-test

data:  par_pre$Tajima and par_post$Tajima
t = 68.637, df = 75560, p-value < 2.2e-16
alternative hypothesis: true mean difference is not equal to 0
95 percent confidence interval:
 0.1683611 0.1782592
sample estimates:
mean difference 
      0.1733102 

mean(par_pre$Tajima)
# 0.5161279
mean(par_post$Tajima)
# 0.3428177