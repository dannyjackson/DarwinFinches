# allele frequency change md 

# filter by linkage
#!/bin/bash

module load python/3.11/3.11.4
module load R
module load plink/1.9
module load samtools

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all

~/programs/angsd/angsd -vcf-gl genolike.bcf -doPlink 2 -out genolike_plink -doGeno -1 -dopost 1 -domaf 1 -doMajorMinor 1 -sites /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs


plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --indep-pairwise 50kb 1 0.5 --out genolike_filtered

plink --tped /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tped --tfam /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_plink.tfam --allow-extra-chr --snps-only 'just-acgt' --extract genolike_filtered.prune.in --out genolike_pruned --make-bed 

# convert to a bed file for use in --sites flag in angsd
awk 'BEGIN{OFS="\t"} {print $1, $4}' genolike_pruned.bim | grep 'NC_' | grep -v 'NC_044601' > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.sites.bed

awk 'NR==FNR { key[$1 FS $2]; next } ($1 FS $2) in key' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/geno_likelihoods/all/genolike_pruned.sites.bed /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs > /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.ldpruned.mafs

~/programs/angsd/angsd sites index /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.ldpruned.mafs

# not working possibly a fault in the sites file
head /xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs

# generate SAF files 

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs

species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=saf_${sp}_pre_pruned \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_pre_pruned.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.4_siteallelefrequency_pruned.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}pre

    sbatch --account=mcnew \
            --job-name=saf_${sp}_post_pruned \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_post_pruned.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.4_siteallelefrequency_pruned.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}post
done






























































































# Testing various scripts for analyzing 

# use SAF files generated during FST analyses
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs

# Estimate Joint 2D SFS
~/programs/angsd/misc/realSFS crapre.saf.idx crapost.saf.idx -P 16 > cra_pre_post.sfs
~/programs/angsd/misc/realSFS forpre.saf.idx forpost.saf.idx > for_pre_post.sfs
~/programs/angsd/misc/realSFS parpre.saf.idx parpost.saf.idx > par_pre_post.sfs


# make mafs files
species=( "cra" "for" "par")
populations=( "pre" "post" )

for sp in "${species[@]}"; do
    for pop in "${populations[@]}"; do
        sbatch --account=mcnew \
                --job-name=posterior.${sp}.${pop} \
                --partition=standard \
                --mail-type=ALL \
                --output=slurm_output/output.posterior.${sp}.${pop}.%j \
                --nodes=1 \
                --ntasks-per-node=4 \
                --time=7:00:00 \
                --mem=100gb \
                ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/deltaAF/PosteriorAF.sh \
                -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
                -s ${sp} -q ${pop}
    done
done

# filter out scaffolds and sex chromosomes

# cra
(zcat cra/cra_pre_postmafs.mafs.gz | head -n 1 && zcat cra/cra_pre_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > cra/cra_pre.mafs.gz
(zcat cra/cra_post_postmafs.mafs.gz | head -n 1 && zcat cra/cra_post_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > cra/cra_post.mafs.gz

# for
(zcat for/for_pre_postmafs.mafs.gz | head -n 1 && zcat for/for_pre_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > for/for_pre.mafs.gz
(zcat for/for_post_postmafs.mafs.gz | head -n 1 && zcat for/for_post_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > for/for_post.mafs.gz

# par
(zcat par/par_pre_postmafs.mafs.gz | head -n 1 && zcat par/par_pre_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > par/par_pre.mafs.gz
(zcat par/par_post_postmafs.mafs.gz | head -n 1 && zcat par/par_post_postmafs.mafs.gz | tail -n +2 | grep 'NC_' | grep -v 'NC_044601') | bgzip > par/par_post.mafs.gz


# compute delta AF and test for significance 
# First, filter snps on linkage
# the linkage disequilibrium between neighboring SNP (using PLINK [98], we removed one SNP from pairs of SNP harboring R2> 0.99 over 20 successive SNP and 50 kb of distance in any sample)


species=( "cra" "for" "par")

for sp in "${species[@]}"; do

    Rscript ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/deltaAF/deltaAF.likelihoodratio.r "/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/${sp}" "${sp}" 

done


CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
cp deltaAF_lrt_significant_q0.20.tsv deltaAF_lrt_significant_q0.20.chrom.tsv
WIN_OUT="deltaAF_lrt_significant_q0.20.chrom.tsv"
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "$WIN_OUT" 
    done < "$CHROM"
fi

CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
cp deltaAF_lrt_significant_q0.08.tsv deltaAF_lrt_significant_q0.08.chrom.tsv
WIN_OUT="deltaAF_lrt_significant_q0.08.chrom.tsv"
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "$WIN_OUT" 
    done < "$CHROM"
fi




CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
cp deltaAF_lrt_significant_p0.001.tsv deltaAF_lrt_significant_p0.001.chrom.tsv
WIN_OUT="deltaAF_lrt_significant_p0.001.chrom.tsv"
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed -i "s/$second/$first/g" "$WIN_OUT" 
    done < "$CHROM"
fi








Rscript ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/deltaAF/signasel_wrapper.r cra "/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/cra" 

species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=signasel.${sp} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.signasel.${sp}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/deltaAF/signasel.sh \
            -s ${sp} -o /xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/${sp}
done


cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af




































# Extract Posterior Allele Frequencies



awk '
BEGIN {
  FS = OFS = "\t"
}

# First file: pre timepoint
FILENAME == ARGV[1] {
  key = $1 "_" $2
  af_pre[key] = $6  # assuming knownEM is in column 6
  next
}

# Second file: post timepoint
FILENAME == ARGV[2] {
  key = $1 "_" $2
  if (key in af_pre) {
    af1 = af_pre[key]
    af2 = $6
    delta = af2 - af1
    print $1, $2, af1, af2, delta
  }
}
' cra_pre.mafs cra_post.mafs > cra_delta_af.txt


CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
WIN_OUT="cra_delta_af.txt"
if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."
    while IFS=',' read -r first second; do
        echo "Replacing $second with $first..."
        sed "s/$second/$first/g" "$WIN_OUT" >> "cra_delta_af.chrom.txt" 
    done < "$CHROM"
fi

echo -e 'chromo\tpos\tpre_af\tpost_af\tneg_log_pvalues_one_tailed' > cra_delta_af.chrom.diff.txt

awk '{ if ($5 > 0) print $0 }' cra_delta_af.chrom.txt >> cra_delta_af.chrom.diff.txt

#!/bin/sh

source ~/programs/DarwinFinches/param_files/cra_params_fst.sh

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/cra

# Run R script for plotting
Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot_bigfile.r" \
    "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "cra_delta_af.chrom.diff.txt" "${WIN}" "delta_af" "cra"

sbatch --account=mcnew \
            --job-name=plot_cra \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.plot_cra.%j \
            --nodes=1 \
            --ntasks-per-node=1 \
            --time=3:00:00 \
            --mem=100gb \
            plot_deltaAF.sh
# to test for sig diff, compute AF per site per indiv, then run t-test?
# what is the distribution of likelihoods genome-wide?


# compute stats 
awk 'NR > 1 {
    x = $5
    sum += x
    sumsq += x * x
    if (x > max) max = x
    n++
} END {
    mean = sum / n
    sd = sqrt((sumsq - sum^2 / n) / (n - 1))
    print "Mean:", mean
    print "SD:", sd
    print "Max:", max
}' cra_delta_af.chrom.diff.txt


# extract top 1%
# Step 1: Extract the 5th column and sort numerically
awk 'NR > 1 { print $5 }' cra_delta_af.chrom.diff.txt | sort -g > values.txt
# Step 2: Find the 99th percentile threshold
threshold=$(awk '{
    a[NR] = $1
}
END {
    n = NR
    cutoff = int(n * 0.99)
    print a[cutoff]
}' values.txt)

# Step 3: Filter out values above the threshold and calculate stats
awk -v thresh="$threshold" 'NR > 1 && $5 <= thresh {
    x = $5
    sum += x
    sumsq += x * x
    if (x > max) max = x
    n++
} END {
    mean = sum / n
    sd = sqrt((sumsq - sum^2 / n) / (n - 1))
    print "Filtered (below 99th percentile)"
    print "Mean:", mean
    print "SD:", sd
    print "Max:", max
}' cra_delta_af.chrom.diff.txt



# compute delta af and test sig
#!/usr/bin/env Rscript

# Load required library
suppressPackageStartupMessages(library(data.table))

# Read command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript compute_delta_af.R crapre.mafs.gz crapost.mafs.gz")
}
pre_file <- "crapre.mafs.gz"
post_file <- "crapost.mafs.gz"

# Read only needed columns: chrom, pos, freq (i.e. minor allele frequency)
cat("Reading pre file...\n")
pre <- fread(cmd = paste("zcat", pre_file), select = c("chromo", "position", "major", "minor", "knownEM", "nInd"))

cat("Reading post file...\n")
post <- fread(cmd = paste("zcat", post_file), select = c("chromo", "position", "major", "minor", "knownEM", "nInd"))

# Merge by chromosome and position
cat("Merging files...\n")
merged <- merge(pre, post, by = c("chromo", "position"))

# Check for mismatches in major and minor alleles across time points
cat("Checking allele identity consistency...\n")
major_mismatch <- merged[major.x != major.y]
minor_mismatch <- merged[minor.x != minor.y]

cat("Major allele mismatches:", nrow(major_mismatch), "\n")
cat("Minor allele mismatches:", nrow(minor_mismatch), "\n")

if (nrow(major_mismatch) > 0 || nrow(minor_mismatch) > 0) {
  cat("Warning: some sites have flipped major/minor alleles.\n")
  cat("Consider resolving these before interpreting delta AF.\n")
} else {
  cat("All major/minor alleles match across timepoints. ✅\n")
}

# Compute ΔAF
cat("Computing delta allele frequency...\n")
merged[, delta_af := knownEM.y - knownEM.x]

# Optional: filter for non-missing
merged <- merged[!is.na(delta_af)]


# test for significance

# Estimate SE using binomial approximation (assuming HWE and diploid)
merged[, pre_se := sqrt(knownEM.x * (1 - knownEM.x) / (2 * nInd.x))]
merged[, post_se := sqrt(knownEM.y * (1 - knownEM.y) / (2 * nInd.y))]
merged[, delta_af := knownEM.y - knownEM.x]
merged[, delta_se := sqrt(pre_se^2 + post_se^2)]

# Z-test
merged[, z := delta_af / delta_se]
merged[, p := 2 * pnorm(-abs(z))]

# Optional: adjust p-values
merged[, q := p.adjust(p, method = "fdr")]

# filter for significant sites
significant_sites <- merged[q < 0.05 & !is.na(q)]

# Output
fwrite(merged, "cra_deltaAF_with_significance.tsv", sep = "\t")
fwrite(significant_sites, "cra_sig_deltaAF.tsv", sep = "\t")

cat("Done.), "\n")




