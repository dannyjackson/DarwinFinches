# script for timesweeper

### must be run on elgato -- breaks with slim > 4
### NOTE: the timesweeper code overwrites the output csv file so you have to edit it to read as follows:
#### file to edit: /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/lib/python3.8/site-packages/timesweeper/find_sweeps_vcf.py
#### predictions.to_csv(outfile, mode="a", header=not os.path.exists(outfile), index=False, sep="\t")
### Compute statistics that are used as parameters in timesweeper
```
# compute average missingness from vcf
module load vcftools
vcftools --vcf /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.vcf --missing-indv --out /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic

awk 'NR > 1 {sum += $5; count++} END {if (count > 0) print sum/count}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.imiss
### 0.0206774

# computed Ne in a different script, see B_PopulationStructureAnalysis/EffectivePopulationSize.md
```
### Load in environment (hashtag is command used to intially create it)
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper

interactive -a mcnew -t 4:00:00 

module load python/3.8/3.8.12
# python3 -m venv --system-site-packages timesweeper_env
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 
```

### cra
#### merge vcf for input
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.unsorted.vcf.gz"
vcf_out_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"

bcftools sort "${vcf_out}" -Oz -o  "${vcf_out_sorted}" 
bcftools index -t "${vcf_out_sorted}" 
# check if output is sorted
bcftools view --no-version "${vcf_out_sorted}"  | grep -v "^#" | sort -k1,1 -k2,2n | diff - <(grep -v "^#" "${vcf_out_sorted}" ) | head

```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out_sorted}" | bgzip > "${vcf_in}"

```
#### Run timesweeper
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

timesweeper sim_custom -y cra.yaml

# remove duplicates from simulated vcfs
cd vcfs
for model in neut ssv; do
  for rep in {0..9}; do
    dir="${model}/${rep}"
    if [ -f "${dir}/merged.vcf" ]; then
      echo "Processing ${dir}/merged.vcf"
      bcftools norm -d both -Ov -o "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
      mv "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
    else
      echo "Warning: ${dir}/merged.vcf not found, skipping."
    fi
  done
done
cd ..

timesweeper condense -o cra.pkl -m 0.02 -y  cra.yaml --hft
timesweeper train -i cra.pkl -y cra.yaml --hft

# first, split vcf into chromosomes
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

# Input VCF and base output directory
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.cleaned.vcf.gz"
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect"
yaml_file="cra.yaml"

# Make output subdirectory for splits
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from VCF header
chroms=$(bcftools view -h "$vcf_in" | grep "^##contig" | sed -E 's/.*ID=([^,]+).*/\1/')

# Loop over chromosomes
for chrom in $chroms; do
    echo "Processing chromosome: $chrom"

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/cra_${chrom}.vcf.gz"

    # Extract chromosome-specific VCF
    bcftools view -r "$chrom" -Oz -o "$vcf_chr" "$vcf_in"

    zgrep -v 'ID=VDB,.*Version=' "${vcf_chr}" | bgzip > "${vcf_chr}.cleaned"

    bcftools index -f "${vcf_chr}.cleaned"

done

rm cra_detect/split_vcfs/cra_NW*

# now run sweepfinder
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

# Input VCF and base output directory
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.cleaned.vcf.gz"
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect"
yaml_file="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra.yaml"

# Make output subdirectory for splits
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from VCF header
CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r chrom; do

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/cra_${chrom}.vcf.gz"

    mkdir -p "${out_base}/chroms/${chrom}/"

    # Run timesweeper detect
    timesweeper detect \
        -i "$vcf_chr" \
        -o "${out_base}/chroms/${chrom}" \
        -y "$yaml_file" \
        --hft

done < "$CHROM"



sbatch --account=mcnew \
    --job-name=detect_cra \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_cra.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_cra.sh
    # Submitted batch job 3947663

echo -e "Chrom\tBP\tPred_Class\tWin_Start\tWin_End\tneut_Prob\tssv_Prob\tssv_selcoeff_pred" > /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_aft.csv

# Append all files excluding their first lines
for file in /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/chroms/*/cra_aft.csv; do
    tail -n +2 "$file" >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_aft.csv
done

cat /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/chroms/*/cra_aft.bed > /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_aft.bed

echo -e "Chrom\tBP\tPred_Class\tWin_Start\tWin_End\tneut_Prob\tssv_Prob\tssv_selcoeff_pred" > /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_hft.csv

# Append all files excluding their first lines
for file in /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/chroms/*/cra_hft.csv; do
    tail -n +2 "$file" >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_hft.csv
done

cat /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/chroms/*/cra_hft.bed > /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_hft.bed

# check if any rows are inconsistent with window orders
awk '$2 < $4 || $4 > $5' cra_aft.csv 
awk '$4 < $2 || $2 > $3' cra_aft.bed
awk '$2 < $4 || $4 > $5' cra_hft.csv 
awk '$4 < $2 || $2 > $3' cra_hft.bed

```

#### Analyze timesweeper results
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/analysis

##### aft
```
# Let's first plot the data to look at the distribution

CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
TS_OUT="cra_aft.csv"
# add new column with numbered chromosomes
awk -F'\t' 'BEGIN {
    FS=OFS="\t"
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[2]] = a[1]
    }
}
NR==1 {
    print "chromo", $0
    next
}
{
    print map[$1], $0
}' "$TS_OUT" > cra_aft.with_chrnum.tsv

module load R

# Load required package
library(ggplot2)

# Set file paths
ts_out <- "cra_aft.with_chrnum.tsv"  # replace with actual path
outdir <- "plots"
dir.create(outdir, showWarnings = FALSE)

# Read in the data
df <- read.table(ts_out, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Histogram of ssv_Prob
p1 <- ggplot(df, aes(x = ssv_Prob)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "black") +
  ggtitle("Histogram of ssv_Prob") +
  xlab("ssv_Prob") + ylab("Count")
ggsave(file.path(outdir, "hist_ssv_prob.png"), p1)

# Histogram of ssv_selcoeff_pred
p2 <- ggplot(df, aes(x = ssv_selcoeff_pred)) +
  geom_histogram(bins = 30, fill = "darkorange", color = "black") +
  ggtitle("Histogram of ssv_selcoeff_pred") +
  xlab("ssv_selcoeff_pred") + ylab("Count")
ggsave(file.path(outdir, "hist_ssv_selcoeff_pred.png"), p2)

# Scatterplot of ssv_Prob vs neut_Prob
p3 <- ggplot(df, aes(x = ssv_Prob, y = neut_Prob)) +
  geom_point(alpha = 0.5) +
  ggtitle("ssv_Prob vs neut_Prob") +
  xlab("ssv_Prob") + ylab("neut_Prob")
ggsave(file.path(outdir, "scatter_ssv_vs_neut_prob.png"), p3)

# Scatterplot of ssv_Prob vs ssv_selcoeff_pred
p4 <- ggplot(df, aes(x = ssv_Prob, y = ssv_selcoeff_pred)) +
  geom_point(alpha = 0.5) +
  ggtitle("ssv_Prob vs ssv_selcoeff_pred") +
  xlab("ssv_Prob") + ylab("ssv_selcoeff_pred")
ggsave(file.path(outdir, "scatter_ssv_vs_ssv_selcoeff_pred.png"), p4)





# analyze aft results
### Identify top 1% of outliers

#### Step 1: Extract the ssv_Prob values (column 7) ignoring the header
awk 'NR > 1 {print $7}' cra_aft.csv | sort -n > ssv_probs.txt

#### Step 2: Find the 99th percentile value
total=$(wc -l < ssv_probs.txt)
index=$(awk -v total="$total" 'BEGIN { printf("%d", total*0.999) }')
cutoff=$(awk -v idx="$index" 'NR == idx {print; exit}' ssv_probs.txt)

echo "Cutoff for top 1% is: $cutoff" # 0.56419164

#### Step 3: Now filter the original file based on this cutoff
awk -v cutoff="$cutoff" 'NR==1 || $7 >= cutoff' cra_aft.csv > cra_aft_outliers.csv

#### make bed file
awk '{print $1,$2-1, $2}' cra_aft_outliers.csv | tail -n +2 | awk '{$1=$1; OFS="\t"}1' > cra_aft_outliers.bed

module load bedtools2 

BED_FILE="cra_aft_outliers.bed"
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
GENES_FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/genes/cra_aft_timesweeper.txt"
GENENAMES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/cra_aft_timesweeper.genenames.txt"
bedtools intersect -a ${GFF} -b ${BED_FILE} -wa > ${GENES_FILE}

grep 'ID\=gene' ${GENES_FILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}

wc -l $GENENAMES
# 2270 

```
##### hft
```
# analyze hft results
module load bedtools2 


CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
TS_OUT="cra_hft.csv"
# add new column with numbered chromosomes
awk -F'\t' 'BEGIN {
    FS=OFS="\t"
    while ((getline < "'$CHROM'") > 0) {
        split($0, a, ",")
        map[a[2]] = a[1]
    }
}
NR==1 {
    print "chromo", $0
    next
}
{
    print map[$1], $0
}' "$TS_OUT" > cra_hft.with_chrnum.tsv

### Identify top 1% of outliers

#### Step 1: Extract the selcoeff values (column 7) ignoring the header
awk 'NR > 1 {print $7}' cra_hft.csv | sort -n > ssv_probs.txt

#### Step 2: Find the 99th percentile value
total=$(wc -l < ssv_probs.txt)
index=$(awk -v total="$total" 'BEGIN { printf("%d", total*0.999) }')
cutoff=$(awk -v idx="$index" 'NR == idx {print; exit}' ssv_probs.txt)

echo "Cutoff for top 1% is: $cutoff" # 0.9607805

#### Step 3: Now filter the original file based on this cutoff
awk -v cutoff="$cutoff" 'NR==1 || $7 >= cutoff' cra_hft.csv > cra_hft_outliers.csv

#### make bed file
awk '{print $1,$2-1, $2}' cra_hft_outliers.csv | tail -n +2 | awk '{$1=$1; OFS="\t"}1' > cra_hft_outliers.bed

BED_FILE="cra_hft_outliers.bed"
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
GENES_FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/genes/cra_hft_timesweeper.txt"
GENENAMES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/cra_hft_timesweeper.genenames.txt"
bedtools intersect -a ${GFF} -b ${BED_FILE} -wa > ${GENES_FILE}
grep 'ID\=gene' ${GENES_FILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${GENENAMES}
wc -l $GENENAMES

# 1710

```



































### for
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

#### merge vcf for input
```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.unsorted.vcf.gz"
vcf_out_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"

bcftools sort "${vcf_out}" -Oz -o  "${vcf_out_sorted}" 
bcftools index -t "${vcf_out_sorted}" 

```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out_sorted}" | bgzip > "${vcf_in}"
```
#### Run timesweeper
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

timesweeper sim_custom -y for.yaml

# remove duplicates from simulated vcfs
cd vcfs
for model in neut ssv; do
  for rep in {0..9}; do
    dir="${model}/${rep}"
    if [ -f "${dir}/merged.vcf" ]; then
      echo "Processing ${dir}/merged.vcf"
      bcftools norm -d both -Ov -o "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
      mv "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
    else
      echo "Warning: ${dir}/merged.vcf not found, skipping."
    fi
  done
done
cd ..

timesweeper condense -o for.pkl -m 0.02 -y  for.yaml --hft
timesweeper train -i for.pkl -y for.yaml --hft


# first, split vcf into chromosomes
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

# Input VCF and base output directory
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.cleaned.vcf.gz"
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for/for_detect"
yaml_file="for.yaml"

bcftools index $vcf_in

# Make output subdirectory for splits
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from reference file
CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r chrom; do
    echo "Processing chromosome: $chrom"

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/for_${chrom}.vcf.gz"

    # Extract chromosome-specific VCF
    bcftools view -r "$chrom" -Oz -o "$vcf_chr" "$vcf_in"

done < "$CHROM"


# now run sweepfinder
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

# Input yaml and base output directory
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for/for_detect"
yaml_file="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for/for.yaml"

# Make output subdirectory for splits
mkdir -p "$out_base"
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from reference file
CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r chrom; do

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/for_${chrom}.vcf.gz"

    mkdir -p "${out_base}/chroms/${chrom}/"

    # Run timesweeper detect
    timesweeper detect \
        -i "$vcf_chr" \
        -o "${out_base}/chroms/${chrom}" \
        -y "$yaml_file" \
        --hft

done < "$CHROM"

sbatch --account=mcnew \
    --job-name=detect_for \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_for.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_for.sh
    # Submitted batch job 3947665
```
### par
#### merge vcf for input
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par
```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.unsorted.vcf.gz"
vcf_out_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"

bcftools sort "${vcf_out}" -Oz -o  "${vcf_out_sorted}" 
bcftools index -t "${vcf_out_sorted}" 

```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out_sorted}" | bgzip > "${vcf_in}"

```
#### Run timesweeper
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par


module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

timesweeper sim_custom -y par.yaml

# remove duplicates from simulated vcfs
cd vcfs
for model in neut ssv; do
  for rep in {0..9}; do
    dir="${model}/${rep}"
    if [ -f "${dir}/merged.vcf" ]; then
      echo "Processing ${dir}/merged.vcf"
      bcftools norm -d both -Ov -o "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
      mv "${dir}/merged.deduped.vcf" "${dir}/merged.vcf"
    else
      echo "Warning: ${dir}/merged.vcf not found, skipping."
    fi
  done
done
cd ..

timesweeper condense -o par.pkl -m 0.02 -y  par.yaml --hft
timesweeper train -i par.pkl -y par.yaml --hft

# first, split vcf into chromosomes
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par

# Input VCF and base output directory
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.cleaned.vcf.gz"
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par/par_detect"
yaml_file="par.yaml"

bcftools index $vcf_in

# Make output subdirectory for splits
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from reference file
CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r chrom; do
    echo "Processing chromosome: $chrom"

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/par_${chrom}.vcf.gz"

    # Extract chromosome-specific VCF
    bcftools view -r "$chrom" -Oz -o "$vcf_chr" "$vcf_in"

    zgrep -v 'ID=VDB,.*Version=' "${vcf_chr}" | bgzip > "${vcf_chr}.cleaned"

    bcftools index -f "${vcf_chr}.cleaned"

done < "$CHROM"


# now run sweepfinder
#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par

# Input yaml and base output directory
out_base="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par/par_detect"
yaml_file="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par/par.yaml"

# Make output subdirectory for splits
mkdir -p "$out_base"
split_dir="${out_base}/split_vcfs"
mkdir -p "$split_dir"

# Get chromosome names from VCF header
CHROM="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/autosomes.txt"

# Loop over chromosomes
while read -r chrom; do

    # Define per-chrom VCF path
    vcf_chr="${split_dir}/par_${chrom}.vcf.gz"

    mkdir -p "${out_base}/chroms/${chrom}/"

    # Run timesweeper detect
    timesweeper detect \
        -i "$vcf_chr" \
        -o "${out_base}/chroms/${chrom}" \
        -y "$yaml_file" \
        --hft

done < "$CHROM"

sbatch --account=mcnew \
    --job-name=detect_par \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_par.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_par.sh
    # Submitted batch job 3947666
```

