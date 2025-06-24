# generate SAF files 
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs

species=( "cra" "for" "par")

species=( "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=saf_${sp}_pre \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_pre.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}pre

    sbatch --account=mcnew \
            --job-name=saf_${sp}_post \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_post.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}post
done

~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh 

/xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs

Set the working directory to fst so that all relevant slurm output files appear here as well:
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst
```
# Run once per species to generate all files

# Define species, environments, and window sizes

# fst
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst

sp=( "cra" )
sp=( "for" )
sp=( "par" )

chmod +x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh 

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
-w 500000 -s 500000

sp=( "par" )

sbatch --account=mcnew \
        --job-name=fst_500000_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.fst_500000_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=3:00:00 \
        --mem=100gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
        -w 500000 -s 500000
        
# Iterate over each combination
species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=fst_500000_${sp} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.fst_500000_${sp}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
            -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
            -w 500000 -s 500000
done

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/allsnps -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa



# Iterate over several window sizes

# Define species, environments, and window sizes
species=( "cra" "for" "par")
window_sizes=( 50000  )


# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        step=$(expr $win / 4)
            ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
            -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
               -w $win -s $step 
    done
done

CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
WIN_OUT="crapre_crapost.50000.fst"
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed "s/$second/$first/g" "$WIN_OUT" >> "${WIN_OUT}.chrom.txt" 
    done < "$CHROM"
fi


window_sizes=( 1 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

        sbatch --account=mcnew \
               --job-name=fst_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.fst_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=1:00:00 \
               --mem=50gb \
                ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
                -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh  \
               -w $win -s $win 
    done
done




sed -i 's/Nsites/Nsites\tfst/' crapre_crapost/50000/crapre_crapost.50000.fst
sed -i 's/Nsites/Nsites\tfst/' forpre_forpost/50000/forpre_forpost.50000.fst
sed -i 's/Nsites/Nsites\tfst/' parpre_parpost/50000/parpre_parpost.50000.fst

# Rerun FST on files with windows that are filtered by depth and mapability
species=( "cra" "for" "par")
# first, add header to input file (double check if necessary before doing so)
for sp in "${species[@]}"; do
    input="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${sp}pre_${sp}post/50000/${sp}pre_${sp}post.50000.fst.depthmapfiltered"

    echo -e "region\tchromo\tposition\tNsites\tfst" | cat - ${input} > tmpfile && mv tmpfile "${input}"
done

# Iterate over each combination
# chmod +x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.filteredfiles.sh
# chmod -x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.filteredfiles.sh

species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    input="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${sp}pre_${sp}post/50000/${sp}pre_${sp}post.50000.fst.depthmapfiltered"

    ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.filteredfiles.sh \
    -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
        -f ${input}
done


# Make violin plots
library(ggplot2)
library(dplyr)

# Read the files
file1 <- read.table("crapre_crapost/50000/crapre_crapost.50000.fst", header = TRUE, stringsAsFactors = FALSE)
file2 <- read.table("forpre_forpost/50000/forpre_forpost.50000.fst", header = TRUE, stringsAsFactors = FALSE)
file3 <- read.table("parpre_parpost/50000/parpre_parpost.50000.fst", header = TRUE, stringsAsFactors = FALSE)

file1 <- read.table("crapre_crapost/50000/crapre_crapost.50000.fst", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(fst >= 0) %>%
  mutate(dataset = "file1")

file2 <- read.table("forpre_forpost/50000/forpre_forpost.50000.fst", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(fst >= 0) %>%
  mutate(dataset = "file2")

file3 <- read.table("parpre_parpost/50000/parpre_parpost.50000.fst", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(fst >= 0) %>%
  mutate(dataset = "file3")

combined <- bind_rows(file1, file2, file3)

combined$dataset <- factor(combined$dataset, 
                           levels = c("file1", "file2", "file3"), 
                           labels = c("cra", "for", "par"))


kruskal_test <- kruskal.test(fst ~ dataset, data = combined)
print(kruskal_test)


        Kruskal-Wallis rank sum test

data:  fst by dataset
Kruskal-Wallis chi-squared = 3675.5, df = 2, p-value < 2.2e-16

# Install FSA
install.packages(dunn.test)
library(dunn.test)

# Post-hoc Dunn's test with Bonferroni correction
dunn_result <- dunn.test(combined$fst, combined$dataset, method = "bonferroni", kw = TRUE)
print(dunn_result)

# 50000 data
data:  fst by dataset
W = 2613407156, p-value < 2.2e-16
alternative hypothesis: true location shift is not equal to 0

# Create and save the plot

p <- ggplot(combined, aes(x = dataset, y = fst, fill = dataset)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  theme_minimal() +
  labs(x = "", y = "FST value", title = "FST Comparison across finch taxa") +
  theme(legend.position = "none")

# Save to file
ggsave("fst_violin_cra_for_par.50000.png", plot = p, width = 9, height = 4, dpi = 300)
