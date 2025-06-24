#!/bin/bash

awk 'NR==1 {next} { 
    sum = 0; 
    for (i = 3; i <= NF; i++) sum += $i; 
    avg = sum / (NF - 2); 
    print $1, $2, avg;
}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_depthstats.txt  > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt 


sbatch --account=mcnew \
--job-name=avgdepthstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.avgdepthstats.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=5:00:00 \
/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/avgdepthstats.sh











#!/bin/bash

# Define window size
WINDOW_SIZE=50000

# Read the input file and process each line
awk -v wsize=$WINDOW_SIZE '
BEGIN { OFS="\t" }
{
    window_start = int(($2 - 1) / wsize) * wsize + 1;
    print $0, $1 "_" window_start;
}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.windowterm.txt


# sites with data as denominator
awk '{sum[$4] += $3; count[$4]++} END {for (w in sum) print w, sum[w]/count[w]}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.windowterm.txt > test_output.txt
# 50000 as denominator
awk '{sum[$4] += $3; count[$4]++} END {for (w in sum) print w, sum[w]/50000}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.windowterm.txt > test_output.2.txt


sed 's/\.1\_/.1 /g' test_output.txt | sort -k1,1 -k2,2n > test_output_sorted.txt
sed 's/\.1\_/.1 /g' test_output.txt | sort -k3,3n | head
sed 's/\.1\_/.1 /g' test_output.txt | sort -k3,3n | tail
awk '{ total += $2 } END { print total/NR }' test_output.txt

R 
# Load data
df <- read.table("test_output_sorted.txt", header=FALSE, sep=" ")

# Assign column names
colnames(df) <- c("chrom", "first_pos", "avg_depth")

# Calculate mean and standard deviation
mean_depth <- mean(df$avg_depth, na.rm=TRUE)
sd_depth <- sd(df$avg_depth, na.rm=TRUE)

# Define cutoffs (±2 SD)
lower_cutoff <- mean_depth - 2 * sd_depth
upper_cutoff <- mean_depth + 2 * sd_depth

# Filter data
filtered_df <- df[df$avg_depth >= lower_cutoff & df$avg_depth <= upper_cutoff, ]

# Save filtered data
write.table(filtered_df, "filtered_windows.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

# Print cutoffs
cat("Lower cutoff:", lower_cutoff, "\n")
cat("Upper cutoff:", upper_cutoff, "\n")

> nrow(df)
[1] 20353
> nrow(filtered_df)
[1] 20230
> lower_cutoff
[1] 3.980108
> upper_cutoff
[1] 12.7701


Some depth updates:

Using unfiltered bam files, my 50kb windows range in average depth from 2.38427 to 241.289, with an average depth of 8.3751. I’m using a cutoff of +/- 2 standard deviations from the mean, which trims 130 windows but keeps 20,230. I think this is really good, definitely will reduce the bias in read depth. 

I haven’t compared it yet to the outlier regions to see if any of the putatively selected regions are from these weirdo regions.

Probably I should also look at mapability in a windowed framework too, right?


# Filter background gene list
awk '{print $9}' /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff | awk -F"[-;]" '{print $2}' | sort -u > genelist.txt


# /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff
# $1 chrom, $4 start, $5 end
awk 'BEGIN {OFS = "\t"} {print $1, $2, $2 + 49999}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/filtered_windows.txt > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/filtered_windows.bed 

# count number of OG genes
cat /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff | grep 'ID\=gene' | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u | wc -l
# 16830

# count number of post-depth-filtering genes
bedtools intersect -a /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/filtered_windows.bed -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff -wb | grep 'ID\=gene' | awk '{OFS = "\t"} {split($12, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u | wc -l
# 15398
# lost 1432 genes due to depth filtering

bedtools intersect -a /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/filtered_windows.bed -b /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.genes.gff -wb | grep 'ID\=gene' | awk '{OFS = "\t"} {split($12, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > genes.windepthfiltered.gff




# Define file paths
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt"

# define bam file
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.bam"

# define output file
# OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.regions.csv"
# this one uses a temp directory to output individual files during threading and then combines later
OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.regions.III.csv"

# define avg depth output file
AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.depth.III.csv"


# run via generalized script on puma
sbatch --account=mcnew \
--job-name=fstavgdepthfastparallel3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fstavgdepthfastparallel3.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows3.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGDEPTH_FILE} -t 94
# 12592202

grep -v "NC_044601" "${AVGDEPTH_FILE}" > "${AVGDEPTH_FILE}.autosomes" 

R 
# Load data
df <- read.table("/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.depth.III.csv.autosomes", header=FALSE, sep="\t")

# Assign column names
colnames(df) <- c("chrom", "win_start", "win_end", "avg_depth", "numsites")

# Calculate mean and standard deviation
mean_depth <- mean(df$avg_depth, na.rm=TRUE)
sd_depth <- sd(df$avg_depth, na.rm=TRUE)

# Define cutoffs (±2 SD)
lower_cutoff <- mean_depth - 2 * sd_depthcd
upper_cutoff <- mean_depth + 2 * sd_depth

# Filter data
filtered_df <- df[df$avg_depth >= lower_cutoff & df$avg_depth <= upper_cutoff, c(1,2,3)]

# Save filtered data
write.table(filtered_df, "50kb.fst.depthfiltered.windows.bam", row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")

# Print cutoffs
cat("Lower cutoff:", lower_cutoff, "\n")
cat("Upper cutoff:", upper_cutoff, "\n")

nrow(df) # 75560
nrow(filtered_df) # 75072
lower_cutoff # 4.43679
upper_cutoff # 12.36991

sort -k1,1 -k2,2n /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kb.fst.depthfiltered.windows.bam > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kb.fst.depthfiltered.windows.sorted.bam

# filter fst file by windows in filtered bam file
BAMFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kb.fst.depthfiltered.windows.sorted.bam"
FSTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.chrom"
OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.depthfiltered"

grep -v "NC_044601" "${FSTFILE}" > "${FSTFILE}.autosomes"

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000
FSTFILE="${FSTFILE}.autosomes"
awk '{print $1, ($2+25000)}' "${BAMFILE}" > "${BAMFILE}.midpos"


awk '
# Load BAMFILE: store all (chrom, midpoint) pairs
FNR==NR {
    bam[$1, $2] = 1
    next
}

# For each line in FSTFILE, keep only if the (chrom, midpoint) pair is in BAMFILE
{
    chrom = $2
    midpoint = $3
    if ((chrom, midpoint) in bam) {
        print
    }
}
' "$BAMFILE" "$FSTFILE" > "$OUTFILE"


wc -l ${OUTFILE}
wc -l ${FSTFILE}




# Replace chromosome names if conversion file is provided
if [ -n "$CHR_FILE" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed "s/$second/$first/g" "${WIN_OUT}.chrom" >> "${WIN_OUT}.chrom.txt" 
    done < "$CHR_FILE"
fi



# z transform windowed data
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" \
    "${OUTDIR}" "${CUTOFF}" "${WIN_OUT}.chrom.txt" "${WIN}" "${POP1}_${POP2}"

Z_OUT="${OUTDIR}/analyses/fst/${POP1}_${POP2}/${POP1}_${POP2}.fst.${WIN}.Ztransformed.csv"

# sed -i 's/\"//g' ${Z_OUT}

# Run R script for plotting
echo "Generating Manhattan plot from ${Z_OUT}..."
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "${WIN}" "fst" "${POP1}" "${POP2}"

echo "Script completed successfully!"












# repeat with mapability scores (see BackgroundGeneList_mapability)
# plot and determine outlier regions
# generate background gene list specific to FST windows









# create with theta windows to confirm script utility
# Define file paths
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt"
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.bam"
OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.regions.csv"
AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.depth.csv"

awk 'BEGIN { OFS="\t" } NR>1 {print $2, ($3-25000), ($3+25000)}' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/50000/crapost.theta.thetasWindow.pestPG > "${WIN_FILE}.numchrom"



CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

# Build the replacement logic into an AWK program
awk -v chr_file="$CHR_FILE" '
BEGIN {
    FS = OFS = "\t"
    while ((getline < chr_file) > 0) {
        split($0, a, ",")
        map[a[1]] = a[2]
    }
}
{
    if ($1 in map) $1 = map[$1]
    print
}
' "${WIN_FILE}.numchrom" | grep 'NC_' > "${WIN_FILE}" 



# run via generalized script on puma
sbatch --account=mcnew \
--job-name=thetasavgdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fstavgdepth.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGDEPTH_FILE} -t 94
# 12603536

























# repeat with RAiSD
CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/chrom_avg_depthstats/chrom_avg_depthstats.txt"
SCRIPT_PATH="/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh"
OUTDIR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd"
BAMSTATS_DIR="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/"
SLURM_OUT="slurm_output"

mkdir -p "${SLURM_OUT}"

# for POP in cra_post cra_pre for_post for_pre par_post par_pre; do

for POP in cra_pre; do
    RAISD_FILE="${OUTDIR}/${POP}/50/RAiSD_Report.${POP}.chromosomes"
    WIN_FILE="${OUTDIR}/${POP}/50/RAiSD_${POP}.windows.bed"
    AVGDEPTH_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.depth.csv"

    awk '{print $1, $3, $4}' "$RAISD_FILE" > "${WIN_FILE}.numchrom"

    awk -v chr_file="$CHR_FILE" '
    BEGIN {
        FS = OFS = " "
        while ((getline < chr_file) > 0) {
            split($0, a, ",")
            map[a[1]] = a[2]
        }
    }
    {
        if ($1 in map) $1 = map[$1]
        print
    }
    ' "${WIN_FILE}.numchrom" | grep 'NC_' | grep -v "NC_044601" > "${WIN_FILE}"

done 

for POP in cra_pre; do
    RAISD_FILE="${OUTDIR}/${POP}/50/RAiSD_Report.${POP}.chromosomes"
    WIN_FILE="${OUTDIR}/${POP}/50/RAiSD_${POP}.windows.bed"
    AVGDEPTH_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.depth.csv"

    sbatch --account=mcnew \
        --job-name=raisdavgdepth_${POP} \
        --partition=standard \
        --mail-type=ALL \
        --output=${SLURM_OUT}/output.raisdavgdepth.${POP}.%j \
        --nodes=1 \
        --ntasks-per-node=94 \
        --time=48:00:00 \
        "${SCRIPT_PATH}" -d "${DEPTH_FILE}" -w "${WIN_FILE}" -a "${AVGDEPTH_FILE}" -t 94
done




# compare win files of fst and tajima to see if they match
wc -l /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.bam
wc -l /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.bam

# are all fst windows the same number? yes! 81185
wc -l crapre_crapost/50000/crapre_crapost.50000.fst.chrom
wc -l forpre_forpost/50000/forpre_forpost.50000.fst.chrom
wc -l parpre_parpost/50000/parpre_parpost.50000.fst.chrom

# are all theta windows the same number? yes! 81680
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/50000/crapre.theta.thetasWindow.pestPG
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/50000/crapost.theta.thetasWindow.pestPG
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/50000/forpre.theta.thetasWindow.pestPG
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/50000/forpost.theta.thetasWindow.pestPG
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/50000/parpre.theta.thetasWindow.pestPG
wc -l /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/50000/parpost.theta.thetasWindow.pestPG

# are all raisd windows the same number? nope. in cra, pre > post. in for and par, pre < post.
wc -l cra_pre/50/RAiSD_Report.cra_pre.chromosomes # 260365
wc -l cra_post/50/RAiSD_Report.cra_post.chromosomes # 166687
wc -l for_pre/50/RAiSD_Report.for_pre.chromosomes # 68055
wc -l for_post/50/RAiSD_Report.for_post.chromosomes # 118331
wc -l par_pre/50/RAiSD_Report.par_pre.chromosomes # 80236
wc -l par_post/50/RAiSD_Report.par_post.chromosomes # 161529