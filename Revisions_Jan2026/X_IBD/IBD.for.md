# Estimating Isolation By Descent
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/IBD/for

module load bcftools
module load samtools
module load htslib

VCF=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.vcf.gz
T1=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/forpresamples.txt
T2=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/forpostsamples.txt
THREADS=16

# filter to IBD snp set

bcftools view -m2 -M2 -v snps -Ou $VCF \
| bcftools filter -e 'F_MISSING>0.05 || MAF<0.05 || QUAL<30' -Oz -o for.ibd.snps.vcf.gz

tabix -p vcf for.ibd.snps.vcf.gz

# check missingness
bcftools view -H for.ibd.snps.vcf.gz | head -n 1
bcftools index -n for.ibd.snps.vcf.gz  # number of records

# impute phasing with beagle
wget https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar

java -Xmx64g -jar beagle.27Feb25.75f.jar \
  gt=for.ibd.snps.vcf.gz \
  out=for.ibd.phased \
  nthreads=16
# output: for.ibd.phased.vcf.gz (largely fully phased)

bcftools query -f'[%GT\t]\n' for.ibd.snps.vcf.gz \
| tr '\t' '\n' \
| awk '
$0 ~ /\|/ {p++}
$0 ~ /\// {u++}
$0 ~ /^\.\.\.$|^\.\.\/\.$|^\.\.\|\.$/ {m++}
END{print "phased(|):",p,"unphased(/):",u,"missing:",m}'

# phased(|): 86902 unphased(/): 534768 missing: 

tabix -p vcf for.ibd.phased.vcf.gz

bcftools query -f'[%GT\t]\n' for.ibd.phased.vcf.gz \
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
bcftools query -f '%CHROM\t%POS\n' for.ibd.phased.vcf.gz \
| awk 'BEGIN{OFS="\t"} { print $1, $1"_"$2, $2/1000000, $2 }' \
> for.ibd.map
```

## Run hapmap
From this github: https://github.com/browning-lab/hap-ibd
```
wget https://faculty.washington.edu/browning/hap-ibd.jar

# no filters
java -Xmx64g -jar hap-ibd.jar \
  gt=for.ibd.phased.vcf.gz \
  map=for.ibd.map \
  out=for.hapibd \
  nthreads=16

```

## Visualize
First, modify files to prepare for plotting:

### Segment length distribution (first sanity plot)
```
zcat for.hapibd.ibd.gz \
| awk '{print $8}' \
| sort -n \
| uniq -c \
| awk '{print $2, $1}' > ibd_lengths.txt
```
### Per-pair total IBD (who shares ancestry?)
```
zcat for.hapibd.ibd.gz \
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
IBD=for.hapibd.ibd.gz
T1=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt 
T2=/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post_pops.txt 

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
zcat for.hapibd.ibd.gz \
 | awk '{print $1,$3}' \
 | sort | uniq -c | sort -nr | head
      5 SRR2917292 SRR2917298
      5 SRR2917290 SRR2917298
      4 SRR2917298 SM1271
      4 SRR2917298 SM1204_all
      4 SRR2917298 SM1083
      4 SM1204_all SM1271
      3 SRR2917298 SM1272
      3 SRR2917298 SM1237
      3 SRR2917293 SRR2917298
      3 SRR2917290 SM1204_all

zcat for.hapibd.ibd.gz |
    grep 'NC_044574'
zcat for.hapibd.ibd.gz \
| awk '{print $5}' \
| sort | uniq -c

      3 NC_044571.1
      6 NC_044572.1
      5 NC_044573.1
     13 NC_044574.1
      1 NC_044576.1
      1 NC_044580.1
     63 NC_044586.1


zcat for.hapibd.ibd.gz \
| awk '$5=="NC_044574.1"' \
| awk 'BEGIN{OFS="\t"} {print $5,$6-1,$7}' \
| sort -k2,2n \
| bedtools merge

NC_044574.1     16543243        29730816
NC_044574.1     37041895        52089042

# In text summary:

# Ask whether long IBD overlaps FST, Tajima's D, composite score
## selection vs demography
I think this looks like a population in decline between T1 and T2, which is driving increased in IBD across all individuals. There is no clear boundary between IBD plots that would show some indv are related while others are not.

# Estimate recent Ne from IBD segment counts
## independent of SFS

# Next steps: drop SM1231 and SM1240 due to high cM values in total IBD with lamich_PL16
