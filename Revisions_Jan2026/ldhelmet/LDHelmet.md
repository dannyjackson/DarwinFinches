# LD Helmet
```
module load vcftools

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/ldhelmet

VCF=/xdisk/mcnew/finches/dannyjackson/finches/analyses/ldhelmet/phased_genomes/cra_phased.vcf.gz

vcftools --gzvcf $VCF --ldhelmet --chr NC_044571.1


```
# Check number of snps: 853178
```
bcftools view -r NC_044571.1 -m2 -M2 -v snps "$VCF" | wc -l
```
# Check missingness: 853665
```
vcftools --gzvcf "$VCF" --chr NC_044571.1 --missing-site --out miss
```
# Are my genotypes phased?
```
bcftools view -r NC_044571.1 "$VCF" | head -n 200 | grep -m 5 -Eo '([01]\|[01]|[01]/[01])'
bcftools view -r NC_044571.1 "$VCF_IND" | head -n 200 | grep -m 5 -Eo '([01]\|[01]|[01]/[01])'
```
# Count number of phased vs unphased sites:
## total 10243980 phased(|) 1074918 unphased(/) 9169062
```
bcftools view -r NC_044571.1 ${VCF} \
  | bcftools query -f'[%GT\n]' \
  | awk '{tot++; if($0~"\\|") ph++; if($0~"/") un++} END{print "total",tot,"phased(|)",ph,"unphased(/)",un}'
```
# Try phasing with BEAGLE instead of whatshap:
## Make conda env:
```
