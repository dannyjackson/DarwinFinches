# RAiSD on imputed phased vcf?
module load micromamba bcftools htslib samtools 
micromamba activate r_puma

mkdir RAiSD
cd RAiSD
wget https://github.com/alachins/raisd/archive/master.zip
unzip master.zip
cd raisd-master
./install-RAiSD.sh

# subset pre and post indv from vcf

bcftools view \
  -S /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapresamples.txt \
  -Oz \
  -o cra.pre.ibd.phased.vcf.gz \
  ../cra.ibd.phased.vcf.gz

gunzip cra.pre.ibd.phased.vcf.gz

bcftools view \
  -S /xdisk/mcnew/finches/dannyjackson/finches/referencelists/crapostsamples.txt \
  -Oz \
  -o cra.post.ibd.phased.vcf.gz \
  ../cra.ibd.phased.vcf.gz

gunzip cra.post.ibd.phased.vcf.gz

REF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna # path to reference genome

~/programs/RAiSD/raisd-master/RAiSD \
    -n cra_pre_raisd \
    -I cra.pre.ibd.phased.vcf \
    -f -O -R -P -C "${REF}" -w 50
    
~/programs/RAiSD/raisd-master/RAiSD \
    -n cra_post_raisd \
    -I cra.post.ibd.phased.vcf \
    -f -O -R -P -C "${REF}" -w 50

# plot in combination
mkdir -p infofiles plots reportfiles

mv RAiSD_Info* infofiles/

mv RAiSD_Plot* plots/

mv RAiSD_Report* reportfiles/

for f in reportfiles/RAiSD_Report.cra_pre_raisd.*; do
    chr=$(basename "$f" | sed 's/.*raisd\.//')
    awk -v chr="$chr" -v OFS="\t" '
    {
        win = int($1 / 50000);
        key = chr FS win;
        sum[key] += $7;
        n[key]   += 1;
    }
    END {
        for (k in sum) {
            split(k,a,FS);
            print a[1], a[2]*50000, (a[2]+1)*50000, sum[k]/n[k];
        }
    }' "$f"
done | sort -k1,1 -k2,2n > cra_pre.mu_50kb.tsv


for f in reportfiles/RAiSD_Report.cra_post_raisd.*; do
    chr=$(basename "$f" | sed 's/.*raisd\.//')
    awk -v chr="$chr" -v OFS="\t" '
    {
        win = int($1 / 50000);
        key = chr FS win;
        sum[key] += $7;
        n[key]   += 1;
    }
    END {
        for (k in sum) {
            split(k,a,FS);
            print a[1], a[2]*50000, (a[2]+1)*50000, sum[k]/n[k];
        }
    }' "$f"
done | sort -k1,1 -k2,2n > cra_post.mu_50kb.tsv


# compute windowed difference 
awk -v OFS="\t" '
# Read pre file first
NR==FNR {
    key = $1 FS $2 FS $3
    pre[key] = $4
    next
}

# Process post file
{
    key = $1 FS $2 FS $3
    if (key in pre) {
        delta = $4 - pre[key]
        print $1, $2, $3, pre[key], $4, delta
    }
}
' cra_pre.mu_50kb.tsv cra_post.mu_50kb.tsv \
| sort -k1,1 -k2,2n > cra_delta_mu_50kb.tsv

# check symmetry of delta
awk '{if ($6 > 0) p++; else n++} END {print p, n}' cra_delta_mu_50kb.tsv



# rename scaffolds for plotting
CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

# chr   start   end   mu_pre   mu_post   delta_mu
FILE="cra_delta_mu_50kb.tsv"
# Read CHROM line by line
while IFS=',' read -r first second; do
    echo "Replacing occurrences of '$second' with '$first' in $FILE"
    sed -i.bak "s/$second/$first/g" "$FILE"
done < "$CHR_FILE"

rm -f "${FILE}.bak"

sed -i 's/ /\t/g' "${FILE}"

