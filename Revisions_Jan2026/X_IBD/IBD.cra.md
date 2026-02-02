# Estimating Isolation By Descent
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/IBD/cra

module load bcftools
module load samtools
module load htslib

VCF=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.vcf.gz
T1=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapresamples.txt
T2=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostsamples.txt
THREADS=16

bcftools view -H "$VCF" | head -n 3

# filter to IBD snp set

bcftools view -m2 -M2 -v snps -Ou $VCF \
| bcftools filter -e 'F_MISSING>0.05 || MAF<0.05 || QUAL<30' -Oz -o cra.ibd.snps.vcf.gz

tabix -p vcf cra.ibd.snps.vcf.gz

# check missingness
bcftools view -H cra.ibd.snps.vcf.gz | head -n 1
bcftools index -n cra.ibd.snps.vcf.gz  # number of records

# impute phasing with beagle
wget https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar

java -Xmx64g -jar beagle.27Feb25.75f.jar \
  gt=cra.ibd.snps.vcf.gz \
  out=cra.ibd.phased \
  nthreads=16
# output: cra.ibd.phased.vcf.gz (largely fully phased)

bcftools query -f'[%GT\t]\n' cra.ibd.snps.vcf.gz \
| tr '\t' '\n' \
| awk '
$0 ~ /\|/ {p++}
$0 ~ /\// {u++}
$0 ~ /^\.\.\.$|^\.\.\/\.$|^\.\.\|\.$/ {m++}
END{print "phased(|):",p,"unphased(/):",u,"missing:",m}'

# phased(|): 86902 unphased(/): 534768 missing: 

tabix -p vcf cra.ibd.phased.vcf.gz

bcftools query -f'[%GT\t]\n' cra.ibd.phased.vcf.gz \
| tr '\t' '\n' \
| awk '
$0 ~ /\|/ {p++}
$0 ~ /\// {u++}
$0 ~ /^\.\.\.$|^\.\.\/\.$|^\.\.\|\.$/ {m++}
END{print "phased(|):",p,"unphased(/):",u,"missing:",m}'
```
# run IBD with HapMap
## Make a simple map (cM = bp/1e6)
```
bcftools query -f '%CHROM\t%POS\n' cra.ibd.phased.vcf.gz \
| awk 'BEGIN{OFS="\t"} { print $1, $1"_"$2, $2/1000000, $2 }' \
> cra.ibd.map
```

## Run hapmap
From this github: https://github.com/browning-lab/hap-ibd
```
wget https://faculty.washington.edu/browning/hap-ibd.jar

# no filters
java -Xmx64g -jar hap-ibd.jar \
  gt=cra.ibd.phased.vcf.gz \
  map=cra.ibd.map \
  out=cra.hapibd \
  nthreads=16


# Conservative
java -Xmx64g -jar hap-ibd.jar \
  gt=cra.ibd.phased.vcf.gz \
  map=cra.ibd.map \
  out=cra.hapibd \
  min-seed=3.0 \
  min-output=3.0 \
  min-markers=200 \
  nthreads=16

# Very recent ancestry only (≈ last ~10–15 generations)
java -Xmx64g -jar hap-ibd.jar \
  gt=cra.ibd.phased.vcf.gz \
  map=cra.ibd.map \
  out=cra.hapibd_5cM \
  min-seed=5.0 \
  min-output=5.0 \
  min-markers=300 \
  nthreads=16
```

## Visualize
First, modify files to prepare for plotting:

### Segment length distribution (first sanity plot)
```
zcat cra.hapibd.ibd.gz \
| awk '{print $8}' \
| sort -n \
| uniq -c \
| awk '{print $2, $1}' > ibd_lengths.txt
```
### Per-pair total IBD (who shares ancestry?)
```
zcat cra.hapibd.ibd.gz \
| awk '{key=$1"__"$3; sum[key]+=$8}
END{for (k in sum) print k, sum[k]}' \
| sort -k2,2nr > ibd_pair_totals.txt
```
### Plot it in R
```
module load micromamba
micromamba activate r_elgato
```

#### Plot 1: IBD segment length distribution → PNG
```
df <- read.table("ibd_lengths.txt")
colnames(df) <- c("cM","count")

png("ibd_length_distribution.png",
    width = 7, height = 5, units = "in", res = 300)

plot(df$cM, df$count, log="y",
     xlab="IBD segment length (cM)",
     ylab="Count (log scale)",
     pch=16)

abline(v=2, col="red", lty=2)

dev.off()
```
#### Plot 2: Per-pair total IBD → PNG
```
df <- read.table("ibd_pair_totals.txt", sep=" ")
colnames(df) <- c("pair","total_cM")

df <- df[order(-df$total_cM),]

png("ibd_pair_totals.png",
    width = 10, height = 5, units = "in", res = 300)

barplot(df$total_cM[1:20],
        names.arg=df$pair[1:20],
        las=2, cex.names=0.6,
        ylab="Total IBD (cM)",
        main="Top IBD-sharing pairs")

dev.off()
```

# Next steps:

# Compare T1–T2 vs T2–T2 IBD totals
## does post-decline sharing increase?

### Compute per-pair total IBD for both classes
```
IBD=cra.hapibd.ibd.gz
T1=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapresamples.txt
T2=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostsamples.txt

# T1–T2
zcat $IBD \
| awk -v T1=$T1 -v T2=$T2 '
BEGIN{
  while ((getline < T1) > 0) pre[$1]=1
  while ((getline < T2) > 0) post[$1]=1
}
{
  if ((pre[$1] && post[$3]) || (pre[$3] && post[$1])) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T1_T2", k, sum[k]
}' > ibd_T1_T2_totals.txt

# T2–T2
zcat $IBD \
| awk -v T2=$T2 '
BEGIN{
  while ((getline < T2) > 0) post[$1]=1
}
{
  if (post[$1] && post[$3] && $1 != $3) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T2_T2", k, sum[k]
}' > ibd_T2_T2_totals.txt

# T1–T1 
zcat $IBD \
| awk -v T1=$T1 '
BEGIN{
  while ((getline < T1) > 0) pre[$1]=1
}
{
  if (pre[$1] && pre[$3] && $1 != $3) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T1_T1", k, sum[k]
}' > ibd_T1_T1_totals.txt

# Combine all for plotting
cat ibd_T1_T1_totals.txt ibd_T1_T2_totals.txt ibd_T2_T2_totals.txt > ibd_pair_totals_by_class.txt

```
Plot:
```
df <- read.table("ibd_pair_totals_by_class.txt")
colnames(df) <- c("class","pair","total_cM")

png("ibd_T1_T2_vs_T2_T2.png", width=6, height=5, units="in", res=300)

boxplot(total_cM ~ class, data=df,
        log="y",
        ylab="Total IBD per pair (cM, log scale)",
        xlab="Pair class",
        main="IBD sharing across vs after decline",
        col=c("grey80","grey50"))

stripchart(total_cM ~ class, data=df,
           vertical=TRUE, method="jitter",
           pch=16, col="black", add=TRUE)

dev.off()

```
### 
zcat cra.hapibd.ibd.gz \
> | awk '{print $1,$3}' \
> | sort | uniq -c | sort -nr | head
      9 lamich_PL16 SM1231
      7 lamich_PL16 SM1240_all
      3 lamich_PL7 SM1266_all
      3 lamich_PL16 lamich_PL4
      2 SM1231 SM1240_all
      2 lamich_PL9 JP9655_all
      2 lamich_PL4 JP4481_all
      2 lamich_PL16 JP4481_all
      2 lamich_PL15 SM1231
      2 lamich_PL15 SM1200_all

zcat cra.hapibd.ibd.gz \
| awk '($1=="lamich_PL16" && $3=="SM1231") || ($1=="SM1231" && $3=="lamich_PL16") {print $5}' \
| sort | uniq -c

      1 NC_044572.1
      8 NC_044577.1

zcat cra.hapibd.ibd.gz \
| awk '($1=="lamich_PL16" && $3=="SM1231") || ($1=="SM1231" && $3=="lamich_PL16") {print $0}' \
| sort | uniq -c | awk '$6=="NC_044572.1"'
#       1 lamich_PL16     1       SM1231  1       NC_044572.1     17135561        19220111        2.084

zcat cra.hapibd.ibd.gz \
| awk '($1=="lamich_PL16" && $3=="SM1231") || ($1=="SM1231" && $3=="lamich_PL16")' \
| awk '$5=="NC_044577.1"' \
| awk 'BEGIN{OFS="\t"} {print $5,$6-1,$7}' \
| sort -k2,2n \
| bedtools merge

NC_044577.1     9376002 12445884
NC_044577.1     13297758        31596510

# In text summary:
Hap-IBD reports IBD at the haplotype level, so overlapping segments with different haplotype combinations (1/2 vs 2/2) represent multiple haplotype–haplotype matches across the same genomic region rather than independent ancestry events.

# Ask whether long IBD overlaps FST, Tajima's D, composite score
## selection vs demography
My personal interpretation of this is that it POTENTIALLY is a haplotype under selection, but it is also potential evidence of recent shared ancestry. If we are to infer selection on a specific gene, multiple broken up haplotypes from different familial lineages with shared SNPs should be changing in frequency over time.

# Estimate recent Ne from IBD segment counts
## independent of SFS

# Next steps: drop SM1231 and SM1240 due to high cM values in total IBD with lamich_PL16

SM1231 and SM1240 (post) have high cM values in total IBD with lamich_PL16 (pre). Can I infer where in the genome these are related?


# Filter and redo
```
bcftools take out samples SM1231 SM1240 from cra.ibd.phased.vcf.gz 

bcftools view \
  -s ^SM1231,SM1240_all \
  -Oz \
  -o cra.ibd.phased.noSM1231_SM1240.vcf.gz \
  cra.ibd.phased.vcf.gz

```
# run IBD with HapMap
## Make a simple map (cM = bp/1e6)
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/IBD/cra/filtered

bcftools query -f '%CHROM\t%POS\n' ../cra.ibd.phased.noSM1231_SM1240.vcf.gz \
| awk 'BEGIN{OFS="\t"} { print $1, $1"_"$2, $2/1000000, $2 }' \
> cra.ibd.map
```

## Run hapmap
From this github: https://github.com/browning-lab/hap-ibd
```
wget https://faculty.washington.edu/browning/hap-ibd.jar

# no filters
java -Xmx64g -jar hap-ibd.jar \
  gt=../cra.ibd.phased.noSM1231_SM1240.vcf.gz \
  map=cra.ibd.map \
  out=cra.hapibd \
  nthreads=16


# Very recent ancestry only (≈ last ~10–15 generations)
java -Xmx64g -jar hap-ibd.jar \
  gt=../cra.ibd.phased.noSM1231_SM1240.vcf.gz \
  map=cra.ibd.map \
  out=cra.hapibd_5cM \
  min-seed=5.0 \
  min-output=5.0 \
  min-markers=300 \
  nthreads=16
```

## Visualize
First, modify files to prepare for plotting:

### Segment length distribution (first sanity plot)
```
zcat cra.hapibd.ibd.gz \
| awk '{print $8}' \
| sort -n \
| uniq -c \
| awk '{print $2, $1}' > ibd_lengths.txt
```
### Per-pair total IBD (who shares ancestry?)
```
zcat cra.hapibd.ibd.gz \
| awk '{key=$1"__"$3; sum[key]+=$8}
END{for (k in sum) print k, sum[k]}' \
| sort -k2,2nr > ibd_pair_totals.txt
```
### Plot it in R
```
module load micromamba
micromamba activate r_elgato
```

#### Plot 1: IBD segment length distribution → PNG
```
df <- read.table("ibd_lengths.txt")
colnames(df) <- c("cM","count")

png("ibd_length_distribution.png",
    width = 7, height = 5, units = "in", res = 300)

plot(df$cM, df$count, log="y",
     xlab="IBD segment length (cM)",
     ylab="Count (log scale)",
     pch=16)

abline(v=2, col="red", lty=2)

dev.off()
```
#### Plot 2: Per-pair total IBD → PNG
```
df <- read.table("ibd_pair_totals.txt", sep=" ")
colnames(df) <- c("pair","total_cM")

df <- df[order(-df$total_cM),]

png("ibd_pair_totals.png",
    width = 10, height = 5, units = "in", res = 300)

barplot(df$total_cM[1:20],
        names.arg=df$pair[1:20],
        las=2, cex.names=0.6,
        ylab="Total IBD (cM)",
        main="Top IBD-sharing pairs")

dev.off()
```

# Next steps:

# Compare T1–T2 vs T2–T2 IBD totals
## does post-decline sharing increase?

### Compute per-pair total IBD for both classes
```
IBD=cra.hapibd.ibd.gz
T1=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapresamples.txt
T2=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostsamples.txt

# T1–T2
zcat $IBD \
| awk -v T1=$T1 -v T2=$T2 '
BEGIN{
  while ((getline < T1) > 0) pre[$1]=1
  while ((getline < T2) > 0) post[$1]=1
}
{
  if ((pre[$1] && post[$3]) || (pre[$3] && post[$1])) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T1_T2", k, sum[k]
}' > ibd_T1_T2_totals.txt

# T2–T2
zcat $IBD \
| awk -v T2=$T2 '
BEGIN{
  while ((getline < T2) > 0) post[$1]=1
}
{
  if (post[$1] && post[$3] && $1 != $3) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T2_T2", k, sum[k]
}' > ibd_T2_T2_totals.txt

# T1–T1 
zcat $IBD \
| awk -v T1=$T1 '
BEGIN{
  while ((getline < T1) > 0) pre[$1]=1
}
{
  if (pre[$1] && pre[$3] && $1 != $3) {
    key=$1"__"$3
    sum[key]+=$8
  }
}
END{
  for (k in sum) print "T1_T1", k, sum[k]
}' > ibd_T1_T1_totals.txt

# Combine all for plotting
cat ibd_T1_T1_totals.txt ibd_T1_T2_totals.txt ibd_T2_T2_totals.txt > ibd_pair_totals_by_class.txt

```
Plot:
```
df <- read.table("ibd_pair_totals_by_class.txt")
colnames(df) <- c("class","pair","total_cM")

png("ibd_T1_T2_vs_T2_T2.png", width=6, height=5, units="in", res=300)

boxplot(total_cM ~ class, data=df,
        log="y",
        ylab="Total IBD per pair (cM, log scale)",
        xlab="Pair class",
        main="IBD sharing across vs after decline",
        col=c("grey80","grey50"))

stripchart(total_cM ~ class, data=df,
           vertical=TRUE, method="jitter",
           pch=16, col="black", add=TRUE)

dev.off()

```
### 
``` 
zcat cra.hapibd.ibd.gz \
 | awk '{print $1,$3}' \
 | sort | uniq -c | sort -nr | head  
```
Output:

    3 lamich_PL7 SM1266_all
    2 lamich_PL9 JP9655_all
    2 lamich_PL4 JP4481_all
    2 lamich_PL16 lamich_PL4
    2 lamich_PL16 JP4481_all
    2 lamich_PL15 SM1200_all
    2 JP5410_all SM1157_all
    1 SM1200_all SM1266_all
    1 lamich_PL9 JP4481_all
    1 lamich_PL7 SM1067_all

``` 
zcat cra.hapibd.ibd.gz \
| awk '{print $5}' \
| sort | uniq -c 
```
Output:

    4 NC_044571.1
    7 NC_044572.1
    6 NC_044573.1
    1 NC_044574.1
    1 NC_044575.1
    1 NC_044576.1
    4 NC_044577.1
    2 NC_044581.1
    1 NC_044586.1
    2 NC_044590.1

### Does TEAD4 occur in one of the haplotypes with high IBD: 
Region: NC_044571.1:87095605- 87147397

```
zcat cra.hapibd.ibd.gz \
| awk '$5=="NC_044571.1"' \
| awk 'BEGIN{OFS="\t"} {print $5,$6-1,$7}' \
| sort -k2,2n \
| bedtools merge
```
Output:

    NC_044571.1     27465805        31334061
    NC_044571.1     47074471        49817902
    NC_044571.1     71834790        76074805 

No, TEAD4 does not.

### Does BMPER occur in one of the haplotypes with high IBD: 
Region: NC_044572.1:46478775-46625955
```
zcat cra.hapibd.ibd.gz \
| awk '$5=="NC_044572.1"' \
| awk 'BEGIN{OFS="\t"} {print $5,$6-1,$7}' \
| sort -k2,2n \
| bedtools merge
```
Output:

    NC_044572.1     6971328 9740480
    NC_044572.1     15032808        17865830
    NC_044572.1     32876664        37927647
    NC_044572.1     71759609        74505693
    NC_044572.1     81547031        84125410
    NC_044572.1     116249860       118345199
    NC_044572.1     122152309       124502570

No, BMPER does not.

### Does NUDT4 occur in one of the haplotypes with high IBD: 
Region: NC_044586.1:28599179-28626256
```
zcat cra.hapibd.ibd.gz \
| awk '$5=="NC_044586.1"' \
| awk 'BEGIN{OFS="\t"} {print $5,$6-1,$7}' \
| sort -k2,2n \
| bedtools merge
```
Output:
  NC_044586.1     10505149        12857525 
No, NUDT4 does not.
