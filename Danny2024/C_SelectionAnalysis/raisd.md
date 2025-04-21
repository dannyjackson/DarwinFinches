# RAiSD.md 

# merge vcfs
# path to lists of samples in each population

sed -i 's/PL/lamich_PL/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre_pops.txt
sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt
sed -i 's/PARV/lamich_PARV/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt
sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt


/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post_pops.txt 


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt
done

for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt
done

for IND in `ls *vcf `
 do
 echo ${IND}
 bgzip ${IND}
 done

# FOR PRE
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

####### HERE #####


sbatch --account=mcnew \
--job-name=indexing \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.indexing.%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=10:00:00 \
/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/indexing.sh 

#!/bin/bash

module load bcftools

# FOR POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done


# CRA PRE
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done


# CRA POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

# PAR PRE

sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

# PAR POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf.gz


#!/bin/bash

# Define arrays for species and time periods
species=("cra" "par" "for")
time_periods=("pre" "post")

# Base directory paths
file_list_base="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering"
output_base="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3"

# Loop through species and time periods
for sp in "${species[@]}"; do
    for tp in "${time_periods[@]}"; do
        input_list="${file_list_base}/${sp}_${tp}.phasedvcfs.txt"
        output_file="${output_base}/${sp}_${tp}.phased.vcf.gz"

        sbatch --account=mcnew \
        --job-name=merge.${sp}.${tp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/merge.${sp}.${tp}.%j \
        --nodes=1 \
        --ntasks-per-node=1 \
        --time=24:00:00 \
        merge.sh -i "$input_list" -o "$output_file"
        
    done
done

#!/bin/bash

while getopts ":i:o:" option; do
    case "${option}" in
        i) INPUT=${OPTARG} ;;
        o) OUTPUT=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

module load bcftools

bcftools merge --file-list "$INPUT" -Oz -o "$OUTPUT"

sbatch --account=mcnew \
    --job-name=raisd.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.raisd.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=10:00:00 \
    ~/programs/raisd_preprocessing.sh 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd

gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf.gz

SPECIES=( "cra_pre" ) #  "cra_post" "for_pre" "for_post" "par_pre" "par_post" 
WINDOW_SIZES=( 50 ) # 10 20 50 100 500 1000 


# Iterate over each combination
for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    sbatch --account=mcnew \
        --job-name=raisd_${WIN}_${SP} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.raisd_${WIN}_${SP}.%j \
        --nodes=1 \
        --time=1:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/raisd/raisd.sh -p ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh -n ${SP} -w ${WIN}
    done 
done

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/raisd/raisd.sh -p ~/programs/DarwinFinches/param_files/cra_pre_params_raisd.sh -n cra_pre -w 50 -r $REF

source /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_base.sh

SPECIES=( "cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post" ) # 
WINDOW_SIZES=( 50 ) # 10 20 50 100 500 1000 


for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    POP=${SP}
    source ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh
    
    cd ${OUTDIR}/analyses/raisd/${SP}/${WIN}

    ~/programs/RAiSD/raisd-master/RAiSD \
        -n "${SP}" \
        -I /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/"${SP}".phased.vcf \
        -f -O -R -P -C "${REF}" -w "${WIN}"
    
    mv *pdf plots
    mv *Report* reportfiles
    mv *Info* infofiles
    done
done


SP="cra_pre"

SPECIES=("cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post") 
# "cra_pre" 
WINDOW_SIZES="50"

# SPECIES=( "for_pre") 

for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh

        POP=${SP}
        echo -e "chrom\tposition\tstart\tend\tVAR\tSFS\tLD\tU" > ${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes


        while read -r SCAFFOLD;
        do 
        awk -v CHROM="$SCAFFOLD" '{print CHROM, $0}' ${OUTDIR}/analyses/raisd/${SP}/${WIN}/reportfiles/RAiSD_Report.${SP}.${SCAFFOLD} >> ${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes

        done < "${OUTDIR}/referencelists/SCAFFOLDS.txt"

        # rename scaffolds for plotting

        echo "Processing CHROM file: $CHR_FILE..."

        FILE="${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes"
        # Read CHROM line by line
        while IFS=',' read -r first second; do
            echo "Replacing occurrences of '$second' with '$first' in $FILE"
            sed -i.bak "s/$second/$first/g" "$FILE"
        done < "$CHR_FILE"

        rm -f "${FILE}.bak"

        sed -i 's/ /\t/g' "${FILE}"


        # z transform U metric


        echo 'z transforming u metric'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" "${OUTDIR}" "${CUTOFF}" "${FILE}" "${WIN}" "${SP}"
        echo 'finished z transforming' 

        WIN_OUT="${OUTDIR}/analyses/raisd/${SP}/${SP}.raisd.${WIN}.Ztransformed.csv"

        sed -i 's/\"//g' $WIN_OUT

        # plot scaffolds
        echo 'visualizing windows'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "raisd" "${SP}"
        echo 'finished windowed plot' 
    done
done

for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh
        
        WIN_OUT="${OUTDIR}/analyses/raisd/${SP}/${SP}.raisd.${WIN}.Ztransformed.csv"

        # plot scaffolds
        echo 'visualizing windows'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "raisd" "${SP}"
        echo 'finished windowed plot' 
    done
done






# plot overlapping


# plot it in R
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd

R 

# CRA

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/finches/dannyjackson/finches"
color2 <- "#4EAFAF"
color1 <- "#082B64"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "cra_post/cra_post.raisd.50.Ztransformed.csv"
input2 <- "cra_pre/cra_pre.raisd.50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "cra_post"
pop2 <- "cra_pre"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)



# Prepare data for plotting
# data 1
cat("Preparing data for plotting...\n")
data1$chromo <- factor(data1$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data1 <- data1 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data1, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data1 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# data 2
cat("Preparing data for plotting...\n")
data2$chromo <- factor(data2$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data2 <- data2 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf2 <- plot_data2 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)

  geom_point(data = plot_data1, aes(x = BPcum, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = plot_data2, aes(x = BPcum, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data1$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "cra.raisd.50.pdf", 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")



# FOR

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/finches/dannyjackson/finches"
color1 <- "#75002B"
color2 <- "#FF817E"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "for_post/for_post.raisd.50.Ztransformed.csv"
input2 <- "for_pre/for_pre.raisd.50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "for_post"
pop2 <- "for_pre"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)



# Prepare data for plotting
# data 1
cat("Preparing data for plotting...\n")
data1$chromo <- factor(data1$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data1 <- data1 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data1, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data1 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# data 2
cat("Preparing data for plotting...\n")
data2$chromo <- factor(data2$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data2 <- data2 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf2 <- plot_data2 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)

  geom_point(data = plot_data1, aes(x = BPcum, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = plot_data2, aes(x = BPcum, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data1$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "for.raisd.50.pdf", 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")



# PAR

# Load required packages, installing if necessary
required_packages <- c("qqman", "hexbin", "readr", "ggrepel", "ggplot2", "dplyr", "RColorBrewer", "data.table")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}
# assign arguments
outdir <- "/xdisk/mcnew/finches/dannyjackson/finches"
color1 <- "#203000"
color2 <- "#A6C965"
cutoff <- as.numeric(0.01)  # Convert to numeric
input1 <- "par_post/par_post.raisd.50.Ztransformed.csv"
input2 <- "par_pre/par_pre.raisd.50.Ztransformed.csv"
win <- "50"
metric <- "raisd"
pop1 <- "par_post"
pop2 <- "par_pre"

# Determine naming convention
pop_name <- ifelse(is.na(pop2), pop1, paste0(pop1, "_", pop2))

# Define parameters
cat("Reading in file...\n")
# Read file
data1 <- fread(input1, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)
data2 <- fread(input2, sep = "\t", na.strings = c("", "NA"), data.table = TRUE)

cat("identifying top snps...\n")
# Identify top SNPs
data_nona <- data1[!is.na(neg_log_pvalues_one_tailed)]
top_snps_count <- round(nrow(data_nona) * cutoff)
cat("identifying top snps 2...\n")
data_nona_sorted <- data_nona %>%
  arrange(desc(neg_log_pvalues_one_tailed)) %>%
  slice_head(n = top_snps_count)

cat("sorting top snps...\n")
# Final sorting
top_snps_dt <- data_nona_sorted[order(chromo, position)]

cat("Get metric cutoff...\n")
# Get metric cutoff
metric_cutoff <- min(top_snps_dt[[metric]], na.rm = TRUE)



# Prepare data for plotting
# data 1
cat("Preparing data for plotting...\n")
data1$chromo <- factor(data1$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data1 <- data1 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data1, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf <- plot_data1 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# data 2
cat("Preparing data for plotting...\n")
data2$chromo <- factor(data2$chromo, levels = c(1, "1A", 2:4, "4A", 5:29, "Z"))

plot_data2 <- data2 %>%
  group_by(chromo) %>%
  summarise(chr_len = max(position)) %>%
  mutate(tot = cumsum(chr_len) - chr_len) %>%
  select(-chr_len) %>%
  left_join(data2, by = "chromo") %>%
  arrange(chromo, position) %>%
  mutate(BPcum = position + tot)

axisdf2 <- plot_data2 %>%
  group_by(chromo) %>%
  summarize(center = mean(BPcum))

# Plot
cat("Generating plot...\n")
ggplot() +
  # Plot data2 (background)

  geom_point(data = plot_data1, aes(x = BPcum, y = !!sym(metric)), 
             color = color1, alpha = 1, size = 0.5) +
  
  # Plot data1 (foreground)
  geom_point(data = plot_data2, aes(x = BPcum, y = !!sym(metric)), 
             color = color2, alpha = 1, size = 0.5) +
  scale_color_manual(values = rep(c(color1, color2), length.out = length(unique(plot_data1$chromo)))) +
  scale_x_continuous(labels = axisdf$chromo, breaks = axisdf$center, guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Chromosome", y = metric) +
  geom_hline(yintercept = metric_cutoff) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


ggsave(filename = "par.raisd.50.pdf", 
       width = 20, height = 5, units = "in")

cat("Script completed successfully!\n")








# Identify outliers and plot using files with windows that are filtered by depth and mapability
species=( "cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post")
# first, add header to input file (double check if necessary before doing so)
for sp in "${species[@]}"; do
    input="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${sp}/50/${sp}.50.raisd.depthmapfiltered"

    echo -e "chrom\tposition\tstart\tend\tVAR\tSFS\tLD\tU" | cat - ${input} > tmpfile && mv tmpfile "${input}"
    rm "${input}.header"
done

# Iterate over each combination

species=("cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post") 
# "cra_pre" 
WINDOW_SIZES="50"

# SPECIES=( "for_pre") 

for WIN in "${WINDOW_SIZES[@]}"; do
    for sp in "${species[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${sp}_params_raisd.sh

        # rename scaffolds for plotting

        echo "Processing CHROM file: $CHR_FILE..."

        FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${sp}/${WIN}/${sp}.${WIN}.raisd.depthmapfiltered"
        cp "${FILE}" "${FILE}.numchrom"

        # Read CHROM line by line
        while IFS=',' read -r first second; do
            sed -i "s/$second/$first/g" "${FILE}.numchrom"
        done < "$CHR_FILE"

        sed -i 's/ /\t/g' "${FILE}"

        # z transform U metric

        echo 'z transforming u metric'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.filteredfiles.r" "${OUTDIR}" "${CUTOFF}" "${FILE}.numchrom"
        echo 'finished z transforming' 

        Z_OUT="${OUTDIR}/analyses/raisd/${sp}/${WIN}/${sp}.${WIN}.raisd.depthmapfiltered.numchrom.Ztransformed.csv"

        sed -i 's/\"//g' $Z_OUT

        # plot scaffolds
        echo 'visualizing windows'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.filteredfiles.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "raisd"
        echo 'finished windowed plot' 
    done
done

species=("cra" "for" "par") 
# "cra_pre" 
WINDOW_SIZES="50"


for WIN in "${WINDOW_SIZES[@]}"; do
    for sp in "${species[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${sp}_pre_params_raisd.sh

        Z_OUT1="${OUTDIR}/analyses/raisd/${sp}_post/${WIN}/${sp}_post.${WIN}.raisd.depthmapfiltered.numchrom.Ztransformed.csv"
        Z_OUT2="${OUTDIR}/analyses/raisd/${sp}_pre/${WIN}/${sp}_pre.${WIN}.raisd.depthmapfiltered.numchrom.Ztransformed.csv"

        OUTFILE="${OUTDIR}/analyses/raisd/${sp}_pre.${sp}_post.${WIN}.raisd.depthmapfiltered.numchrom.Ztransformed.overlapping.png"

        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.filteredfiles.overlapping.r" "${OUTDIR}" "${COLOR2}" "${COLOR1}" "${CUTOFF}" "${Z_OUT1}" "${Z_OUT2}" "raisd" "${OUTFILE}"
   done
done

