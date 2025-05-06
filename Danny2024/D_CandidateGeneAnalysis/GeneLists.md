# running interactively because it is a quick script
chmod +x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/
### FST
```
# Define species, environments, and window sizes
# Identify genes in windows
species=( "cra" "for" "par")
window_sizes=( 50000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
    -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${sp}pre_${sp}post/50000/${sp}pre_${sp}post.50000.fst.depthmapfiltered.numchrom.Ztransformed.csv.outlier.csv \
    -n ${sp} \
    -m fst \
    -w ${win}

    done
done
```
#### Make background gene list
```
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file

MAPFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.fst.map.filtered.bam"

DEPTHFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.fst.depth.filtered.bam"

MAPBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/fst.map.windows.bed"
DEPTHBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/fst.depth.windows.bed"

BACKGROUNDFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/genes/background_genes.fst.txt"

BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.fst.txt"


tr ' ' '\t' < ${MAPFILE} > ${MAPBED}
tr ' ' '\t' < ${DEPTHFILE} > ${DEPTHBED}

bedtools intersect -a ${GFF} -b ${MAPBED} > "${BACKGROUNDFILE}.justmap"
bedtools intersect -a "${BACKGROUNDFILE}.justmap" -b ${DEPTHBED} > ${BACKGROUNDFILE}


grep 'ID\=gene' ${BACKGROUNDFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${BACKGROUNDGENES}


grep 'ID\=gene' ${GFF} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/allgenes.txt

```
#### Analyze overlapping genes (on local computer)
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GeneLists

R

library(ggvenn)

df_cra<-read.csv("cra.fst.50000kb.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("for.fst.50000kb.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("par.fst.50000kb.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))
write(intersect_all, file="intersection.fst.all.csv")

pdf(file = "../plots/intersection.fst.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

#### Clean up panther output (run files through panther and save locally in github repository)
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GOterms/pantheroutput

##### cra
tail -n +12 cra.fst.50000kb.panther.txt > cra.fst.50000kb.panther.clean.txt

sed -i '' 's/cra.fst.50000kb.genenames.txt //g' cra.fst.50000kb.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.fst.50000kb.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.fst.50000kb.panther.clean

head -1 cra.fst.50000kb.panther.clean.txt > cra.fst.50000kb.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' cra.fst.50000kb.panther.clean.txt >> cra.fst.50000kb.panther.clean.p01.txt

##### for
tail -n +12 for.fst.50000kb.panther.txt > for.fst.50000kb.panther.clean.txt

sed -i '' 's/for.fst.50000kb.genenames.txt //g' for.fst.50000kb.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' for.fst.50000kb.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r for.fst.50000kb.panther.clean

head -1 for.fst.50000kb.panther.clean.txt > for.fst.50000kb.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' for.fst.50000kb.panther.clean.txt >> for.fst.50000kb.panther.clean.p01.txt

##### par
```
tail -n +12 par.fst.50000kb.panther.txt > par.fst.50000kb.panther.clean.txt

sed -i '' 's/par.fst.50000kb.genenames.txt //g' par.fst.50000kb.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.fst.50000kb.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.fst.50000kb.panther.clean

head -1 par.fst.50000kb.panther.clean.txt > par.fst.50000kb.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' par.fst.50000kb.panther.clean.txt >> par.fst.50000kb.panther.clean.p01.txt

# grep -o 'GO:[0-9]\+' cra.fst.50000kb.panther.clean.txt | head > go_terms.txt
```

### Tajima diff
```
~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
-p ~/programs/DarwinFinches/param_files/crapost_params_tajima.sh \
-i /xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/cra.Tajima_50000.bottom.outlier.csv \
-n cra_bottom \
-m tajima_diff \
-w 50000
mv gene_names/crapost.tajima_diff.50000kb.genenames.txt gene_names/crapost.tajima_diff.50000kb.bottom.genenames.txt

~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
-p ~/programs/DarwinFinches/param_files/crapost_params_tajima.sh \
-i /xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/cra.Tajima_50000.top.outlier.csv \
-n cra_top \
-m tajima_diff \
-w 50000
mv gene_names/crapost.tajima_diff.50000kb.genenames.txt gene_names/crapost.tajima_diff.50000kb.top.genenames.txt
```
#### chromosome 1

grep 'NC_044571.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr1.bed


GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr1.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr1.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr1.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr1.genenames.txt

#### chromosome 9
grep 'NC_044579.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr9.bed

GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr9.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr9.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr9.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr9.genenames.txt

#### chromosome 11
grep 'NC_044581.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr11.bed

GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr11.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr11.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr11.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr11.genenames.txt

#### post depth and map filtering
species=( "cra" "for" "par")

for sp in "${species[@]}"; do

  ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
  -p ~/programs/DarwinFinches/param_files/${sp}post_params_tajima.sh \
  -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}_tajimadiff.depthmapfiltered.txt.numchrom.Ztransformed.csv.bottom.outlier.csv \
  -n ${sp}_bottom \
  -m tajima_diff \
  -w 50000

  mv /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/${sp}post.tajima_diff.50000kb.genenames.txt /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/${sp}.tajima_diff.50000kb.bottom.genenames.txt

  ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
  -p ~/programs/DarwinFinches/param_files/${sp}post_params_tajima.sh \
  -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}_tajimadiff.depthmapfiltered.txt.numchrom.Ztransformed.csv.top.outlier.csv \
  -n ${sp}_top \
  -m tajima_diff \
  -w 50000

  mv /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/${sp}post.tajima_diff.50000kb.genenames.txt /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/${sp}.tajima_diff.50000kb.top.genenames.txt

done 

#### Make background gene list
```
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file

MAPFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.thetas.map.filtered.bam"

DEPTHFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.thetas.depth.filtered.bam"

MAPBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/thetas.map.windows.bed"
DEPTHBED="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/bed/thetas.depth.windows.bed"

BACKGROUNDFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/genes/background_genes.thetas.txt"

BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.thetas.txt"


tr ' ' '\t' < ${MAPFILE} > ${MAPBED}
tr ' ' '\t' < ${DEPTHFILE} > ${DEPTHBED}

bedtools intersect -a ${GFF} -b ${MAPBED} > "${BACKGROUNDFILE}.justmap"
bedtools intersect -a "${BACKGROUNDFILE}.justmap" -b ${DEPTHBED} > ${BACKGROUNDFILE}


grep 'ID\=gene' ${BACKGROUNDFILE} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ${BACKGROUNDGENES}


grep 'ID\=gene' ${GFF} | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/allgenes.txt

```

#### Clean up panther output (run files through panther and save locally in github repository)
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GOterms/pantheroutput

##### cra
species=( "cra" "for" "par")

for sp in "${species[@]}"; do

  ###### bottom
  tail -n +12 ${sp}.tajima_diff.50000kb.bottom.panther.txt > ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  sed -i '' "s/${sp}.tajima_diff.50000kb.bottom.genenames.txt //g" ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  sed -i '' "s/(raw P-value)/P-value/g" ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r ${sp}.tajima_diff.50000kb.bottom.panther.clean

  head -1 ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt > ${sp}.tajima_diff.50000kb.bottom.panther.clean.p01.txt

  awk 'BEGIN {FS = "\t"} $7 < 0.01' ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt >> ${sp}.tajima_diff.50000kb.bottom.panther.clean.p01.txt

  ###### top

  tail -n +12 ${sp}.tajima_diff.50000kb.top.panther.txt > ${sp}.tajima_diff.50000kb.top.panther.clean.txt

  sed -i '' "s/${sp}.tajima_diff.50000kb.top.genenames.txt //g" ${sp}.tajima_diff.50000kb.top.panther.clean.txt

  sed -i '' "s/(raw P-value)/P-value/g" ${sp}.tajima_diff.50000kb.top.panther.clean.txt

  Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r ${sp}.tajima_diff.50000kb.top.panther.clean

  head -1 ${sp}.tajima_diff.50000kb.top.panther.clean.txt > ${sp}.tajima_diff.50000kb.top.panther.clean.p01.txt

  awk 'BEGIN {FS = "\t"} $7 < 0.01' ${sp}.tajima_diff.50000kb.top.panther.clean.txt >> ${sp}.tajima_diff.50000kb.top.panther.clean.p01.txt

 
done


# raisd


species=( "cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post" )
window_sizes=( 50 )



## Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    echo ${sp}
    ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/DarwinFinches/param_files/${sp}_params_raisd.sh \
    -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${sp}/${win}/${sp}.${win}.raisd.depthmapfiltered.numchrom.Ztransformed.csv.outlier.csv \
    -n ${sp} \
    -m raisd \
    -w ${win}

    done
done


chmod +x ~/programs/DarwinFinches/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh
species=( "cra" "for" "par" )
window_sizes=( 50 )

for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/DarwinFinches/Genomics-Main/D_GeneVisualization/D2_uniquegenes.sh \
    -p ~/programs/DarwinFinches/param_files/params_base.sh \
    -i ${sp}_post \
    -q ${sp}_pre \
    -m raisd \
    -w ${win}

    done
done

## make background gene lists

species=( "cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post" )


for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=raisd_background_${sp} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.raisd_background_${sp}.%j \
            --nodes=1 \
            --ntasks-per-node=1 \
            --time=12:00:00 \
            --mem=30 \
            ~/programs/DarwinFinches/Genomics-Main/D_GeneVisualization/D2_raisd_backgroundgenes.sh \
            -p ${sp}
done


## filter background gene lists locally
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GeneLists

library(ggvenn)

df_cra_pre<-read.csv("background_genes.cra_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_cra_post<-read.csv("background_genes.cra_post.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_for_pre<-read.csv("background_genes.for_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for_post<-read.csv("background_genes.for_post.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

df_par_pre<-read.csv("background_genes.par_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par_post<-read.csv("background_genes.par_post.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]


cra_intersect <- intersect(df_cra_pre, df_cra_post)
write(cra_intersect, file="background_genes.cra_all.raisd.txt")

for_intersect <- intersect(df_for_pre, df_for_post)
write(for_intersect, file="background_genes.for_all.raisd.txt")

par_intersect <- intersect(df_par_pre, df_par_post)
write(par_intersect, file="background_genes.par_all.raisd.txt")


## for unique genes in cra post that were not in cra pre
## 1. Genes in df_cra_post that are not in df_cra_pre
df_cra_genes<-read.csv("cra_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_cra_background<-read.csv("background_genes.cra_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_cra_post <- intersect(df_cra_genes, df_cra_background)

write(unique_cra_post, file="cra.raisd.txt")

## 2. Genes in df_for_post that are not in df_for_pre
df_for_genes<-read.csv("for_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for_background<-read.csv("background_genes.for_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_for_post <- intersect(df_for_genes, df_for_background)

write(unique_for_post, file="for.raisd.txt")


## 3. Genes in df_par_post that are not in df_par_pre
df_par_genes<-read.csv("par_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par_background<-read.csv("background_genes.par_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_par_post <- intersect(df_par_genes, df_par_background)

write(unique_par_post, file="par.raisd.txt")

## evaluate intersection of all raisd genes

cra_for_raisd_intersect <- intersect(unique_cra_post, unique_for_post)
cra_par_raisd_intersect <- intersect(unique_cra_post, unique_par_post)
par_for_raisd_intersect <- intersect(unique_par_post, unique_for_post)

## intersection of all background genes
two_raisd_intersect_background <- intersect(cra_intersect, for_intersect)
all_raisd_intersect_background <- intersect(two_raisd_intersect_background, par_intersect)
write(all_raisd_intersect_background, file="background.all.raisd.txt")

## intersection of all candidate genes
all_raisd_intersect <- intersect(cra_for_raisd_intersect, unique_par_post)
filtered_par_intersect <- intersect(all_raisd_intersect, all_raisd_intersect_background)
write(all_raisd_intersect, file="all.raisd.txt")

## visualize overlapping genes

gene_sets <- list(
  cra= unique_cra_post,
  "for" = unique_for_post,
  par = unique_par_post
)


pdf(file = "intersection.all.raisd.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

##### cra
tail -n +12 cra.raisd.panther.txt > cra.raisd.panther.clean.txt

sed -i '' 's/cra.raisd.txt //g' cra.raisd.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.raisd.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.raisd.panther.clean

head -1 cra.raisd.panther.clean.txt > cra.raisd.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' cra.raisd.panther.clean.txt >> cra.raisd.panther.clean.p01.txt

##### for
tail -n +12 for.raisd.panther.txt > for.raisd.panther.clean.txt

sed -i '' 's/for.raisd.txt //g' for.raisd.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' for.raisd.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r for.raisd.panther.clean

head -1 for.raisd.panther.clean.txt > for.raisd.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' for.raisd.panther.clean.txt >> for.raisd.panther.clean.p01.txt

##### par
tail -n +12 par.raisd.panther.txt > par.raisd.panther.clean.txt

sed -i '' 's/par.raisd.txt //g' par.raisd.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.raisd.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.raisd.panther.clean

head -1 par.raisd.panther.clean.txt > par.raisd.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' par.raisd.panther.clean.txt >> par.raisd.panther.clean.p01.txt


##### all: no genes in overlap

# Delta AF
## Windowed

```
#### Analyze overlapping genes (on local computer)
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GeneLists

R

library(ggvenn)

df_cra<-read.csv("cra.delta_af.windows.1percent.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("for.delta_af.windows.1percent.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("par.delta_af.windows.1percent.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))
write(intersect_all, file="intersection.deltaAF.windowed.all.csv")

# Combine all genes and tabulate their frequency
all_genes <- c(df_cra, df_for, df_par)
gene_counts <- table(all_genes)

# Keep genes present in at least 2 of the 3 lists
genes_in_two_or_more <- names(gene_counts[gene_counts >= 2])

# Write to file
write(genes_in_two_or_more, file = "intersection.deltaAF.windowed.any2.csv")

pdf(file = "../plots/intersection.deltaAF.windowed.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()
```
#### Clean up panther output (run files through panther and save locally in github repository)
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GOterms/pantheroutput

##### cra
tail -n +12 cra.delta_af.windows.1percent.panther.txt > cra.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/cra.delta_af.windows.1percent.genenames.txt //g' cra.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.delta_af.windows.1percent.panther.clean.txt


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.delta_af.windows.1percent.panther.clean

head -1 cra.delta_af.windows.1percent.panther.clean.txt > cra.delta_af.windows.1percent.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' cra.delta_af.windows.1percent.panther.clean.txt >> cra.delta_af.windows.1percent.panther.clean.p01.txt

##### for

tail -n +12 for.delta_af.windows.1percent.panther.txt > for.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/for.delta_af.windows.1percent.genenames.txt //g' for.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' for.delta_af.windows.1percent.panther.clean.txt


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r for.delta_af.windows.1percent.panther.clean

head -1 for.delta_af.windows.1percent.panther.clean.txt > for.delta_af.windows.1percent.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' for.delta_af.windows.1percent.panther.clean.txt >> for.delta_af.windows.1percent.panther.clean.p01.txt


##### par

tail -n +12 par.delta_af.windows.1percent.panther.txt > par.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/par.delta_af.windows.1percent.genenames.txt //g' par.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.delta_af.windows.1percent.panther.clean.txt


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.delta_af.windows.1percent.panther.clean

head -1 par.delta_af.windows.1percent.panther.clean.txt > par.delta_af.windows.1percent.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' par.delta_af.windows.1percent.panther.clean.txt >> par.delta_af.windows.1percent.panther.clean.p01.txt


##### intersction
tail -n +12 intersectionany2.delta_af.windows.1percent.panther.txt > intersectionany2.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/intersection.deltaAF.windowed.any2.csv //g' intersectionany2.delta_af.windows.1percent.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' intersectionany2.delta_af.windows.1percent.panther.clean.txt


Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r intersectionany2.delta_af.windows.1percent.panther.clean

head -1 intersectionany2.delta_af.windows.1percent.panther.clean.txt > intersectionany2.delta_af.windows.1percent.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' intersectionany2.delta_af.windows.1percent.panther.clean.txt >> intersectionany2.delta_af.windows.1percent.panther.clean.p01.txt





# Make plots of FDR passing genes
## RAISD
##### Only CRA had any significant GO terms
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.raisd.panther.clean.fdr.txt > cra.raisd.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.raisd.panther.clean.fdr.GO

## FST
### cra.fst.50000kb.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.fst.50000kb.panther.clean.fdr.txt > cra.fst.50000kb.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fst.50000kb.panther.clean.fdr.GO

### for.fst.50000kb.panther.clean.fdr.txt -- only has one GO term
grep -o 'GO:[0-9]\{7\}' pantheroutput/for.fst.50000kb.panther.clean.fdr.txt > for.fst.50000kb.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.fst.50000kb.panther.clean.fdr.GO


## Tajima
### cra.tajima_diff.50000kb.bottom.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.tajima_diff.50000kb.bottom.panther.clean.fdr.txt > cra.tajima_diff.50000kb.bottom.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.tajima_diff.50000kb.bottom.panther.clean.fdr.GO

### cra.tajima_diff.50000kb.top.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.tajima_diff.50000kb.top.panther.clean.fdr.txt > cra.tajima_diff.50000kb.top.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.tajima_diff.50000kb.top.panther.clean.fdr.GO

### for.tajima_diff.50000kb.top.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/for.tajima_diff.50000kb.top.panther.clean.fdr.txt > for.tajima_diff.50000kb.top.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.tajima_diff.50000kb.top.panther.clean.fdr.GO





## generate all gene lists 
cat cra* | sort -u > cra.fstraisd.genenames.txt
cat for* | sort -u > for.fstraisd.genenames.txt
cat par* | sort -u > par.fstraisd.genenames.txt

## prepare GO term lists for visualization
```
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/cra.fst.50000kb.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/cra_post.raisd.50kb.unique.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/cra.fstraisd.GO

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/for.fst.50000kb.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/for_post.raisd.50kb.unique.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/for.fstraisd.GO

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/par.fst.50000kb.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/par_post.raisd.50kb.unique.genenames
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/par.fstraisd.GO

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/intersection.all
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/intersection.cra_for
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/intersection.cra_par
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/intersection.for_par
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/intersection.for_par_cra

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/cra.raisd.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/for.raisd.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/par.raisd.GO
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r pantheroutput/all.raisd.GO

grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.fst.50000kb.genenames.fdr.txt > cra.fst.50000kb.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra_post.raisd.50kb.unique.genenames.fdr.txt > cra_post.raisd.50kb.unique.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.fstraisd.GO.fdr.txt > cra.fstraisd.GO.txt

grep -o 'GO:[0-9]\{7\}' pantheroutput/for.fst.50000kb.genenames.fdr.txt > for.fst.50000kb.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/for_post.raisd.50kb.unique.genenames.fdr.txt > for_post.raisd.50kb.unique.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/for.fstraisd.GO.fdr.txt > for.fstraisd.GO.txt

grep -o 'GO:[0-9]\{7\}' pantheroutput/par.fst.50000kb.genenames.fdr.txt > par.fst.50000kb.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/par_post.raisd.50kb.unique.genenames.fdr.txt > par_post.raisd.50kb.unique.genenames.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/par.fstraisd.GO.fdr.txt > par.fstraisd.GO.txt

grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.all.fdr.txt > intersection.all.GO.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.cra_for.fdr.txt > intersection.cra_for.GO.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.cra_par.fdr.txt > intersection.cra_par.GO.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.for_par.fdr.txt > intersection.for_par.GO.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/intersection.for_par_cra.fdr.txt > intersection.for_par_cra.GO.txt


## visualize
## Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fst.50000kb.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra_post.raisd.50kb.unique.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fstraisd.GO

# Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.fst.50000kb.genenames.fdr
# Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for_post.raisd.50kb.unique.genenames.fdr

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r par.fst.50000kb.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r par_post.raisd.50kb.unique.genenames.fdr
```















# Analyze overlap between methods
library(ggvenn)

df_cra<-read.csv("GeneLists/cra.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("GeneLists/for.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("GeneLists/par.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

## compare across species
gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))
write(intersect_all, file="GeneLists/intersection.all.csv")

intersect_cra_for = intersect(df_cra, df_for)
write(intersect_cra_for, file="GeneLists/intersection.cra_for.csv")

intersect_cra_par = intersect(df_cra, df_par)
write(intersect_cra_par, file="GeneLists/intersection.cra_par.csv")

intersect_for_par = intersect(df_for, df_par)
write(intersect_for_par, file="GeneLists/intersection.for_par.csv")


combined_unique <- unique(c(intersect_cra_for, intersect_cra_par, intersect_for_par))
write(combined_unique, file="GeneLists/intersection.for_par_cra.csv")

pdf(file = "intersection.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()



## test for significance in overlap
### cra vs for
background<-read.csv("GeneLists/backgroundgenelist.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
N=length(background) # 16628
k = length(intersect_cra_for) # 73
m = length(df_cra) # 590
n = length(df_for) # 466

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value < 2.2e-16
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 4.250684 7.351356
sample estimates:
odds ratio 
  5.619609 

### cra vs par
background<-read.csv("GeneLists/backgroundgenelist.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
N=length(background) # 16628
k = length(intersect_cra_par) # 53
m = length(df_cra) # 590
n = length(df_par) # 668

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 4.268e-08
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 1.809303 3.328471
sample estimates:
odds ratio 
  2.474883 

### for vs par
background<-read.csv("GeneLists/backgroundgenelist.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
N=length(background) # 16628
k = length(intersect_for_par) # 52
m = length(df_for) # 466
n = length(df_par) # 668

fisher.test(matrix(c(k, m-k, n-k, N-m-n+k), nrow = 2))


        Fisher's Exact Test for Count Data

data:  matrix(c(k, m - k, n - k, N - m - n + k), nrow = 2)
p-value = 2.408e-11
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 2.302470 4.289474
sample estimates:
odds ratio 
  3.169484 




## first, make universal background gene dataset
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/

## Set directory and file pattern
FILES=/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background*

## Take the first file as the starting point
cp $(ls $FILES | head -n 1) temp_intersect.txt

## Loop through all remaining files and intersect
for f in $(ls $FILES | tail -n +2); do
    echo "Intersecting with $f..."
    grep -Fxf temp_intersect.txt "$f" > temp_new.txt
    mv temp_new.txt temp_intersect.txt
done

## Final result
mv temp_intersect.txt backgroundgenes_allfiles.txt



# PCA of all stats across genes
## CRA
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_analyses/cra
### define files
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/backgroundgenes_allfiles.txt"
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/cradiff.txt"
RAISD="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/cra_post/50/cra_post.50.raisd.depthmapfiltered"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/cra/deltaAF_lrt_full.windowedavg.tsv"
TIMESWEEPER="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect/cra_aft.csv"



### Output file
OUTFILE="cra.genes_with_stats.tsv"

echo -e "gene\tchr\tstart\tend\tavg_fst\tavg_tajima\tavg_deltaaf\tavg_raisd" > "$OUTFILE"

### Read each gene
while read gene; do
    # Find gene info
    gff_info=$(grep -F "ID=gene-$gene" "$GFF" | grep 'ID=gene' | head -n 1)
    
    # If gene not found, skip
    if [ -z "$gff_info" ]; then
        echo "Warning: $gene not found in GFF."
        continue
    fi

    # Extract chromosome, start, end
    chr=$(echo "$gff_info" | awk '{print $1}')
    start=$(echo "$gff_info" | awk '{print $4}')
    end=$(echo "$gff_info" | awk '{print $5}')

    # Overlap if: window_end >= gene_start && window_start <= gene_end

    # FST (assumes columns: chrom, win_start, win_end, ..., stat_col)
    fst_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
    BEGIN { FS=OFS="\t" }
    NR > 1 {
        # extract actual start and end from the third parenthetical group
        match($1, /\([0-9]+,[0-9]+\)\([0-9]+,[0-9]+\)\(([0-9]+),([0-9]+)\)/, arr)
        win_start = arr[1]
        win_end = arr[2]
        if ($2 == chr && win_end >= gstart && win_start <= gend) {
            print $5
        }
    }
' "$FST")

    # Tajima (assumes: chrom, position (mid), TajimaD)
    tajima_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr {
            win_start = $2 - 25000; win_end = $2 + 25000;
            if (win_end >= gstart && win_start <= gend) print $3;
        }
    ' "$TAJIMA")

    # deltaAF (assumes: chrom, win_start, win_end, delta_af)
    deltaaf_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $3 >= gstart && $2 <= gend {print $4}
    ' "$DELTAAF")

    # RAiSD (assumes: chrom, position, win_start, win_end, ..., ..., ..., U)
    raisd_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $4 >= gstart && $3 <= gend {print $8}
    ' "$RAISD")

    # Averages
    avg_fst=$(echo "$fst_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_tajima=$(echo "$tajima_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_deltaaf=$(echo "$deltaaf_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_raisd=$(echo "$raisd_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')

    # Output line
    echo -e "${gene}\t${chr}\t${start}\t${end}\t${avg_fst}\t${avg_tajima}\t${avg_deltaaf}\t${avg_raisd}" >> "$OUTFILE"

done < "$BACKGROUNDGENES"



### filter out NA values and unidentified (LOC) genes

grep -v 'NA' cra.genes_with_stats.tsv | grep -v 'LOC115' > cra.genes_with_stats.filtered.tsv

### add column with relevant GO term status (TRUE or FALSE is related to angiogenesis)

awk 'NR==FNR {genes[$1]; next} $1 == "gene" {print $0, "\tangiogenesis_status"; next} $1 in genes {$0 = $0 "\tTRUE"} !($1 in genes) {$0 = $0 "\tFALSE"} {print}' angiogenesis.txt cra.genes_with_stats.filtered.tsv > cra.genes_with_stats.filtered.angiogenesis.tsv


```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("cra.genes_with_stats.filtered.angiogenesis.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("avg_fst", "avg_tajima", "avg_deltaaf", "avg_raisd")]), ]
df_clean$avg_deltaaf <- abs(df_clean$avg_deltaaf)

### Select only the numeric columns
stat_matrix <- df_clean[, c("avg_fst", "avg_tajima", "avg_deltaaf", "avg_raisd")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


            avg_fst avg_tajima avg_deltaaf avg_raisd
avg_fst       1.000      0.152       0.307     0.084
avg_tajima    0.152      1.000       0.120     0.110
avg_deltaaf   0.307      0.120       1.000    -0.077
avg_raisd     0.084      0.110      -0.077     1.000



### PCA on fst, tajima, and deltaaf
pca_input <- df_clean[, c("avg_fst", "avg_tajima", "avg_deltaaf")]
pca_result <- prcomp(pca_input, scale. = TRUE)

### Scores = PCA coordinates
scores <- as.data.frame(pca_result$x)  # same number of rows as df_clean
scores$gene <- df_clean$gene  # keep gene names
scores$high_PC1 <- scores$PC1 >= quantile(scores$PC1, 0.99)  # this now matches correctly

### Loadings = variable contribution directions
loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale





### Add PC scores back to df
df_clean$PC1 <- pca_result$x[, 1]
df_clean$PC2 <- pca_result$x[, 2]

df_clean$PC1_2 <- df_clean$PC1*df_clean$PC2

### Identify top 1% by PC1
cutoff <- quantile(df_clean$PC1, probs = 0.99)
df_clean$high_PC1 <- df_clean$PC1 >= cutoff

### Identify top 1% by PC2 (low)
cutoff <- quantile(df_clean$PC2, probs = 0.01)
df_clean$high_PC2 <- df_clean$PC2 <= cutoff

### Identify top 1% by PC1_2 (high)
cutoff <- quantile(df_clean$PC1_2, probs = 0.01)
df_clean$high_PC1_2 <- df_clean$PC1_2 <= cutoff

### Write top genes to file
# top_genes_pc1 <- df_clean[df_clean$high_PC1, ]
# top_genes_pc2 <- top_genes_pc1[top_genes_pc1$high_PC2, ]
top_genes_pc1_2 <- df_clean[df_clean$high_PC1_2, ]$gene

write.table(top_genes_pc1_2, file = "highPC1_lowPC2_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$high_PC1] <- "PC1"
df_clean$highlight_group[df_clean$high_PC2] <- "PC2"
df_clean$highlight_group[df_clean$high_PC1_2] <- "PC1_2"
df_clean$highlight_group[df_clean$high_PC1 & df_clean$high_PC2] <- "PC1 & PC2"

# Convert to factor for custom color ordering
df_clean$highlight_group <- factor(df_clean$highlight_group, levels = c("None", "PC1", "PC2", "PC1_2", "PC1 & PC2"))

pdf(file = "genes.pca.pdf", width = 8, height = 6, useDingbats = FALSE)

ggplot(df_clean, aes(x = PC1, y = PC2)) +
  # First, plot all points with a default color
  geom_point(aes(color = "Other"), size = 1, alpha = 0.8) +
  
  # Now, plot angiogenesis points (red) on top
  geom_point(data = subset(df_clean, angiogenesis_status == TRUE), 
             aes(color = "angiogenesis"), size = 1, alpha = 1) +
  
  scale_color_manual(values = c("Other" = "gray80", 
                                "angiogenesis" = "red")) +  # Set red color for angiogenesis points


  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", size = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +

  theme_minimal(base_size = 14) +
  labs(title = "PCA of Gene Stats with Loadings and Top 1% Genes on PC1 and PC2",
      x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)"),
      color = "Angiogenesis Status")

dev.off()



tail -n +12 cra.PCA.panther.fdr.txt > cra.PCA.panther.fdr.clean.txt

sed -i '' 's/highPC1_lowPC2_genes.tsv //g' cra.PCA.panther.fdr.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.PCA.panther.fdr.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.PCA.panther.fdr.clean

head -1 cra.PCA.panther.fdr.clean.txt > cra.PCA.panther.fdr.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' cra.PCA.panther.fdr.clean.txt >> cra.PCA.panther.fdr.clean.p01.txt


















### Read and clean input

df_clean$selection_score <- with(df_clean,
                                       avg_fst * avg_deltaaf * (-avg_tajima))

### Rank descending (higher score = stronger selection signal)
df_clean <- df_clean[order(-df_clean$selection_score), ]


### Identify top 1% by PC1
cutoff <- quantile(df_clean$selection_score, probs = 0.99)
df_clean$high_selectionscore <- df_clean$selection_score >= cutoff

### Write top genes to file
selected_genes <- df_clean[df_clean$high_selectionscore, ]

### Write ranked results
write.table(selected_genes,
            "genes_ranked_by_selection_score.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

awk '{print $1}' genes_ranked_by_selection_score.tsv | tail -n +2 > genes_ranked_by_selection_score.genenames.tsv


CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
WIN_OUT=$TAJIMA

if [ -n "$CHROM" ]; then
    echo "Replacing chromosome names based on conversion file..."

    TMP_OUT="${WIN_OUT}.tmp"

    awk -v convfile="$CHROM" '
    BEGIN {
        FS = OFS = " "  # Space-separated input and output
        while ((getline < convfile) > 0) {
            split($0, a, ",")
            map[a[1]] = a[2]
        }
        close(convfile)
    }
    {
        if ($1 in map) {
            $1 = map[$1]
        }
        print
    }
    ' "$WIN_OUT" > "$TMP_OUT"

    mv "$TMP_OUT" "$WIN_OUT"
fi


```
### playing around with composite statistic
```
#### Load required package
library(dplyr)
library(ggplot2)
#### Read the data
gene_stats <- read.table("cra.genes_with_stats.noLOC.tsv", header = TRUE, sep = "\t")
gene_stats <- na.omit(gene_stats)
#### Create standardized (z-score) columns
gene_stats <- gene_stats %>%
  mutate(
    z_deltaaf = scale(avg_deltaaf)[,1],            # Strong directional change
    z_tajima = scale(-avg_tajima)[,1],             # More negative D = stronger sweep-like signal
    z_fst = scale(avg_fst)[,1]                     # Temporal differentiation
  )

# Compute composite selection score
# Weighted: deltaAF (0.4), Tajima's D (0.4), FST (0.2)
gene_stats <- gene_stats %>%
  mutate(
    composite_score = 0.4 * z_deltaaf + 0.4 * z_tajima + 0.2 * z_fst
  )

# Identify top 1% genes
threshold <- quantile(gene_stats$composite_score, 0.99)
gene_stats <- gene_stats %>%
  mutate(top_candidate = composite_score >= threshold)

# Optional: rank genes by score (descending)
gene_stats <- gene_stats %>%
  arrange(desc(composite_score))

# Save to file
write.table(gene_stats, "composite_selection_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



# Plot 1: Histogram of composite scores
pdf(file = "composite_hist.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = composite_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Composite Selection Score",
       x = "Composite Score", y = "Number of Genes") +
  theme_minimal()
dev.off()

# Plot 2: Composite Score vs Standardized ΔAF, label top genes
pdf(file = "composite_vs_deltaAF.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = z_deltaaf, y = composite_score)) +
  geom_point(aes(color = top_candidate), alpha = 0.7) +
  geom_text(data = subset(gene_stats, top_candidate), aes(label = gene),
            vjust = -0.5, size = 3) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  labs(title = "Composite Score vs Standardized ΔAF",
       x = "Z(ΔAF)", y = "Composite Score") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()

gene_stats$abs_deltaaf <- abs(gene_stats$avg_deltaaf)

# Plot 3: Standardized ΔAF, label top genes
pdf(file = "tajima_vs_deltaAF.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = abs_deltaaf, y = avg_tajima)) +
  geom_point(aes(color = top_candidate), alpha = 0.2) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  labs(title = "Tajima's D vs Absolute Delta AF",
       x = "Abs Delta AF", y = "Tajima") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()


gene_stats %>%
  arrange(desc(abs_deltaaf)) %>%
  slice(1:5)
```

## FOR

### define files
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/backgroundgenes_allfiles.txt"
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/forpre_forpost/50000/forpre_forpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/fordiff.txt"
RAISD="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/for_post/50/for_post.50.raisd.depthmapfiltered"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/for/deltaAF_lrt_full.windowedavg.tsv"
TIMESWEEPER="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for/for_detect/for_aft.csv"
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

### Convert numeric chromosomes in Tajima file using CHROM mapping
TMP_TAJIMA="tajima.tmp"

awk 'FNR==NR {map[$1]=$2; next}
     BEGIN {OFS="\t"}
     FNR == 1 {print; next} 
     ($1 in map) { $1 = map[$1]; print }' FS="," "$CHROM" FS="[ \t]+" "$TAJIMA" >> tajima.tmp

TAJIMA=$TMP_TAJIMA

### Output file
OUTFILE="for.genes_with_stats.tsv"
echo -e "gene\tchr\tstart\tend\tavg_fst\tavg_tajima\tavg_deltaaf\tavg_raisd" > "$OUTFILE"

### Read each gene
while read gene; do
    # Find gene info
    gff_info=$(grep -F "ID=gene-$gene" "$GFF" | grep 'ID=gene' | head -n 1)
    
    # If gene not found, skip
    if [ -z "$gff_info" ]; then
        echo "Warning: $gene not found in GFF."
        continue
    fi

    # Extract chromosome, start, end
    chr=$(echo "$gff_info" | awk '{print $1}')
    start=$(echo "$gff_info" | awk '{print $4}')
    end=$(echo "$gff_info" | awk '{print $5}')

    # Overlap if: window_end >= gene_start && window_start <= gene_end

    # FST (assumes columns: chrom, win_start, win_end, ..., stat_col)
    fst_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
    BEGIN { FS=OFS="\t" }
    NR > 1 {
        # extract actual start and end from the third parenthetical group
        match($1, /\([0-9]+,[0-9]+\)\([0-9]+,[0-9]+\)\(([0-9]+),([0-9]+)\)/, arr)
        win_start = arr[1]
        win_end = arr[2]
        if ($2 == chr && win_end >= gstart && win_start <= gend) {
            print $5
        }
    }
' "$FST")

    # Tajima (assumes: chrom, position (mid), TajimaD)
    tajima_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr {
            win_start = $2 - 25000; win_end = $2 + 25000;
            if (win_end >= gstart && win_start <= gend) print $3;
        }
    ' "$TAJIMA")

    # deltaAF (assumes: chrom, win_start, win_end, delta_af)
    deltaaf_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $3 >= gstart && $2 <= gend {print $4}
    ' "$DELTAAF")

    # RAiSD (assumes: chrom, position, win_start, win_end, ..., ..., ..., U)
    raisd_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $4 >= gstart && $3 <= gend {print $8}
    ' "$RAISD")

    # Averages
    avg_fst=$(echo "$fst_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_tajima=$(echo "$tajima_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_deltaaf=$(echo "$deltaaf_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_raisd=$(echo "$raisd_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')

    # Output line
    echo -e "${gene}\t${chr}\t${start}\t${end}\t${avg_fst}\t${avg_tajima}\t${avg_deltaaf}\t${avg_raisd}" >> "$OUTFILE"

done < "$BACKGROUNDGENES"


### filter out NA values and unidentified (LOC) genes

grep -v 'NA' for.genes_with_stats.tsv | grep -v 'LOC115' > for.genes_with_stats.filtered.tsv

### add column with relevant GO term status (TRUE or FALSE is related to angiogenesis)

awk 'NR==FNR {genes[$1]; next} $1 == "gene" {print $0, "\tangiogenesis_status"; next} $1 in genes {$0 = $0 "\tTRUE"} !($1 in genes) {$0 = $0 "\tFALSE"} {print}' angiogenesis.txt for.genes_with_stats.filtered.tsv > for.genes_with_stats.filtered.angiogenesis.tsv


```
### Load libraries
library(ggplot2)

### Read and clean input
df <- read.delim("for.genes_with_stats.filtered.angiogenesis.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("avg_fst", "avg_tajima", "avg_deltaaf", "avg_raisd")]), ]
df_clean$avg_deltaaf <- abs(df_clean$avg_deltaaf)

### Select only the numeric columns
stat_matrix <- df_clean[, c("avg_fst", "avg_tajima", "avg_deltaaf", "avg_raisd")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))


            avg_fst avg_tajima avg_deltaaf avg_raisd
avg_fst       1.000      0.027       0.432     0.054
avg_tajima    0.027      1.000      -0.018    -0.065
avg_deltaaf   0.432     -0.018       1.000     0.188
avg_raisd     0.054     -0.065       0.188     1.000



### PCA on fst, tajima, and deltaaf
pca_input <- df_clean[, c("avg_fst", "avg_tajima", "avg_deltaaf")]
pca_result <- prcomp(pca_input, scale. = TRUE)

### Scores = PCA coordinates
scores <- as.data.frame(pca_result$x)  # same number of rows as df_clean
scores$gene <- df_clean$gene  # keep gene names
scores$high_PC1 <- scores$PC1 >= quantile(scores$PC1, 0.99)  # this now matches correctly

### Loadings = variable contribution directions
loadings <- as.data.frame(pca_result$rotation)
loadings$variable <- rownames(loadings)
arrow_scale <- 3
loadings$xend <- loadings$PC1 * arrow_scale
loadings$yend <- loadings$PC2 * arrow_scale





### Add PC scores back to df
df_clean$PC1 <- pca_result$x[, 1]
df_clean$PC2 <- pca_result$x[, 2]

df_clean$PC1_2 <- df_clean$PC1*df_clean$PC2

### Identify top 1% by PC1
cutoff <- quantile(df_clean$PC1, probs = 0.99)
df_clean$high_PC1 <- df_clean$PC1 >= cutoff

### Identify top 1% by PC2 (low)
cutoff <- quantile(df_clean$PC2, probs = 0.01)
df_clean$high_PC2 <- df_clean$PC2 <= cutoff

### Identify top 1% by PC1_2 (high)
cutoff <- quantile(df_clean$PC1_2, probs = 0.01)
df_clean$high_PC1_2 <- df_clean$PC1_2 <= cutoff

### Write top genes to file
# top_genes_pc1 <- df_clean[df_clean$high_PC1, ]
# top_genes_pc2 <- top_genes_pc1[top_genes_pc1$high_PC2, ]
top_genes_pc1_2 <- df_clean[df_clean$high_PC1_2, ]$gene

write.table(top_genes_pc1_2, file = "highPC1_lowPC2_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


df_clean$highlight_group <- "None"
df_clean$highlight_group[df_clean$high_PC1] <- "PC1"
df_clean$highlight_group[df_clean$high_PC2] <- "PC2"
df_clean$highlight_group[df_clean$high_PC1_2] <- "PC1_2"
df_clean$highlight_group[df_clean$high_PC1 & df_clean$high_PC2] <- "PC1 & PC2"

# Convert to factor for custom color ordering
df_clean$highlight_group <- factor(df_clean$highlight_group, levels = c("None", "PC1", "PC2", "PC1_2", "PC1 & PC2"))

pdf(file = "genes.pca.pdf", width = 8, height = 6, useDingbats = FALSE)

ggplot(df_clean, aes(x = PC1, y = PC2)) +
  # First, plot all points with a default color
  geom_point(aes(color = "Other"), size = 1, alpha = 0.8) +
  
  # Now, plot angiogenesis points (red) on top
  geom_point(data = subset(df_clean, angiogenesis_status == TRUE), 
             aes(color = "angiogenesis"), size = 1, alpha = 1) +
  
  scale_color_manual(values = c("Other" = "gray80", 
                                "angiogenesis" = "red")) +  # Set red color for angiogenesis points


  # Arrows for variable loadings
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", size = 1) +

  # Text labels for loading vectors
  geom_text(data = loadings,
            aes(x = xend, y = yend, label = variable),
            color = "black", size = 4, vjust = -0.5, inherit.aes = FALSE) +

  theme_minimal(base_size = 14) +
  labs(title = "PCA of Gene Stats with Loadings and Top 1% Genes on PC1 and PC2",
      x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)"),
      color = "Angiogenesis Status")

dev.off()



### Ensure complete cases and abs(deltaAF)
df_clean <- df[complete.cases(df[, c("avg_fst", "avg_tajima", "avg_deltaaf")]), ]
df_clean$avg_deltaaf <- abs(df_clean$avg_deltaaf)

### Select only the numeric columns
stat_matrix <- df_clean[, c("avg_fst", "avg_tajima", "avg_deltaaf")]

### Compute Pearson correlation matrix
cor_matrix <- cor(stat_matrix, method = "pearson")
print(round(cor_matrix, 3))



gene    avg_fst   avg_tajima  avg_deltaaf
OTULIN  0.131336  -1.29574    0.0616755
TCFL5   0.08312   -0.709907   0.0721721
ITGA2B  0.284453  -0.546081   0.0174979
LETM2   0.098342  -0.042659   0.0608066

### max fst:        0.339061
### min tajima:     -3.61341
### max tajima:     4.32922
### max delta af:   0.104388


### Read and clean input
df <- read.delim("cra.genes_with_stats.noLOC.tsv", header = TRUE, stringsAsFactors = FALSE)
df_clean <- df[complete.cases(df[, c("avg_fst", "avg_tajima", "avg_deltaaf")]), ]
df_clean$avg_deltaaf <- abs(df_clean$avg_deltaaf)

df_clean$selection_score <- with(df_clean,
                                       avg_fst * avg_deltaaf * (-avg_tajima))

### Rank descending (higher score = stronger selection signal)
df_clean <- df_clean[order(-df_clean$selection_score), ]


### Identify top 1% by PC1
cutoff <- quantile(df_clean$selection_score, probs = 0.99)
df_clean$high_selectionscore <- df_clean$selection_score >= cutoff

### Write top genes to file
selected_genes <- df_clean[df_clean$high_selectionscore, ]

### Write ranked results
write.table(selected_genes,
            "genes_ranked_by_selection_score.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

awk '{print $1}' genes_ranked_by_selection_score.tsv | tail -n +2 > genes_ranked_by_selection_score.genenames.tsv




```
### playing around with composite statistic
```
#### Load required package
library(dplyr)
library(ggplot2)
#### Read the data
gene_stats <- read.table("cra.genes_with_stats.noLOC.tsv", header = TRUE, sep = "\t")
gene_stats <- na.omit(gene_stats)
#### Create standardized (z-score) columns
gene_stats <- gene_stats %>%
  mutate(
    z_deltaaf = scale(avg_deltaaf)[,1],            # Strong directional change
    z_tajima = scale(-avg_tajima)[,1],             # More negative D = stronger sweep-like signal
    z_fst = scale(avg_fst)[,1]                     # Temporal differentiation
  )

# Compute composite selection score
# Weighted: deltaAF (0.4), Tajima's D (0.4), FST (0.2)
gene_stats <- gene_stats %>%
  mutate(
    composite_score = 0.4 * z_deltaaf + 0.4 * z_tajima + 0.2 * z_fst
  )

# Identify top 1% genes
threshold <- quantile(gene_stats$composite_score, 0.99)
gene_stats <- gene_stats %>%
  mutate(top_candidate = composite_score >= threshold)

# Optional: rank genes by score (descending)
gene_stats <- gene_stats %>%
  arrange(desc(composite_score))

# Save to file
write.table(gene_stats, "composite_selection_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



# Plot 1: Histogram of composite scores
pdf(file = "composite_hist.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = composite_score)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Composite Selection Score",
       x = "Composite Score", y = "Number of Genes") +
  theme_minimal()
dev.off()

# Plot 2: Composite Score vs Standardized ΔAF, label top genes
pdf(file = "composite_vs_deltaAF.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = z_deltaaf, y = composite_score)) +
  geom_point(aes(color = top_candidate), alpha = 0.7) +
  geom_text(data = subset(gene_stats, top_candidate), aes(label = gene),
            vjust = -0.5, size = 3) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  labs(title = "Composite Score vs Standardized ΔAF",
       x = "Z(ΔAF)", y = "Composite Score") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()

gene_stats$abs_deltaaf <- abs(gene_stats$avg_deltaaf)

# Plot 3: Standardized ΔAF, label top genes
pdf(file = "tajima_vs_deltaAF.pdf", width = 6, height = 6, useDingbats=FALSE)

ggplot(gene_stats, aes(x = abs_deltaaf, y = avg_tajima)) +
  geom_point(aes(color = top_candidate), alpha = 0.2) +
  scale_color_manual(values = c("gray70", "firebrick")) +
  labs(title = "Tajima's D vs Absolute Delta AF",
       x = "Abs Delta AF", y = "Tajima") +
  theme_minimal() +
  theme(legend.position = "none")

dev.off()


gene_stats %>%
  arrange(desc(abs_deltaaf)) %>%
  slice(1:5)
```

## PAR


### define files
GFF="/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff" # path to gff file
BACKGROUNDGENES="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/backgroundgenes_allfiles.txt"
FST="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/parpre_parpost/50000/parpre_parpost.50000.fst.autosomes"
TAJIMA="/xdisk/mcnew/finches/dannyjackson/finches/analyses/Tajima/pardiff.txt"
RAISD="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/par_post/50/par_post.50.raisd.depthmapfiltered"
DELTAAF="/xdisk/mcnew/finches/dannyjackson/finches/analyses/delta_af/par/deltaAF_lrt_full.windowedavg.tsv"
TIMESWEEPER="/xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par/par_detect/par_aft.csv"
CHROM="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"

### Convert numeric chromosomes in Tajima file using CHROM mapping
TMP_TAJIMA="tajima.tmp"

awk 'FNR==NR {map[$1]=$2; next}
     BEGIN {OFS="\t"}
     FNR == 1 {print; next} 
     ($1 in map) { $1 = map[$1]; print }' FS="," "$CHROM" FS="[ \t]+" "$TAJIMA" >> tajima.tmp

TAJIMA=$TMP_TAJIMA

### Output file
OUTFILE="par.genes_with_stats.2.tsv"
echo -e "gene\tchr\tstart\tend\tavg_fst\tavg_tajima\tavg_deltaaf\tavg_raisd" > "$OUTFILE"

### Read each gene
while read gene; do
    # Find gene info
    gff_info=$(grep -F "ID=gene-$gene" "$GFF" | grep 'ID=gene' | head -n 1)
    
    # If gene not found, skip
    if [ -z "$gff_info" ]; then
        echo "Warning: $gene not found in GFF."
        continue
    fi

    # Extract chromosome, start, end
    chr=$(echo "$gff_info" | awk '{print $1}')
    start=$(echo "$gff_info" | awk '{print $4}')
    end=$(echo "$gff_info" | awk '{print $5}')

    # Overlap if: window_end >= gene_start && window_start <= gene_end

    # FST (assumes columns: chrom, win_start, win_end, ..., stat_col)
    fst_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
    BEGIN { FS=OFS="\t" }
    NR > 1 {
        # extract actual start and end from the third parenthetical group
        match($1, /\([0-9]+,[0-9]+\)\([0-9]+,[0-9]+\)\(([0-9]+),([0-9]+)\)/, arr)
        win_start = arr[1]
        win_end = arr[2]
        if ($2 == chr && win_end >= gstart && win_start <= gend) {
            print $5
        }
    }
' "$FST")

    # Tajima (assumes: chrom, position (mid), TajimaD)
    tajima_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr {
            win_start = $2 - 25000; win_end = $2 + 25000;
            if (win_end >= gstart && win_start <= gend) print $3;
        }
    ' "$TAJIMA")

    # deltaAF (assumes: chrom, win_start, win_end, delta_af)
    deltaaf_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $3 >= gstart && $2 <= gend {print $4}
    ' "$DELTAAF")

    # RAiSD (assumes: chrom, position, win_start, win_end, ..., ..., ..., U)
    raisd_vals=$(awk -v chr="$chr" -v gstart="$start" -v gend="$end" '
        BEGIN{OFS="\t"}
        $1 == chr && $4 >= gstart && $3 <= gend {print $8}
    ' "$RAISD")

    # Averages
    avg_fst=$(echo "$fst_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_tajima=$(echo "$tajima_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_deltaaf=$(echo "$deltaaf_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')
    avg_raisd=$(echo "$raisd_vals" | awk '{sum+=$1} END {if (NR>0) print sum/NR; else print "NA"}')

    # Output line
    echo -e "${gene}\t${chr}\t${start}\t${end}\t${avg_fst}\t${avg_tajima}\t${avg_deltaaf}\t${avg_raisd}" >> "$OUTFILE"

done < "$BACKGROUNDGENES"




## Background gene list

FST_BACKGROUND="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.fst.txt"
THETA_BACKGROUND="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.thetas.txt"
PCA_BACKGROUND="/xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.pca.txt"

sort /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.fst.txt > fst.sorted.txt
sort /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/gene_names/background_genes.thetas.txt > thetas.sorted.txt

comm -12 $FST_BACKGROUND $THETA_BACKGROUND > $PCA_BACKGROUND

### it's keeping only the fst genes, so I can just use that as my background gene list

## Panther analysis
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GOterms/PCA/pantheroutput

### CRA 
#### high
tail -n +12 cra.high.PCA.panther.txt > cra.high.PCA.panther.clean.txt

sed -i '' 's/cra.high_PC1_2.genenames.txt //g' cra.high.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.high.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.high.PCA.panther.clean

#### highest
tail -n +12 cra.highest.PCA.panther.txt > cra.highest.PCA.panther.clean.txt

sed -i '' 's/cra.highest_PC1_2.genenames.txt //g' cra.highest.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' cra.highest.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r cra.highest.PCA.panther.clean

### FOR

#### high
tail -n +12 for.high.PCA.panther.txt > for.high.PCA.panther.clean.txt

sed -i '' 's/for.high_PC1_2.genenames.txt //g' for.high.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' for.high.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r for.high.PCA.panther.clean

#### highest
tail -n +12 for.highest.PCA.panther.txt > for.highest.PCA.panther.clean.txt

sed -i '' 's/for.highest_PC1_2.genenames.txt //g' for.highest.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' for.highest.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r for.highest.PCA.panther.clean

### PAR

#### high
tail -n +12 par.high.PCA.panther.txt > par.high.PCA.panther.clean.txt

sed -i '' 's/par.high_PC1_2.genenames.txt //g' par.high.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.high.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.high.PCA.panther.clean

#### highest
tail -n +12 par.highest.PCA.panther.txt > par.highest.PCA.panther.clean.txt

sed -i '' 's/par.highest_PC1_2.genenames.txt //g' par.highest.PCA.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.highest.PCA.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.highest.PCA.panther.clean

## Overlapping gene analysis
cd /Users/danjack/Documents/Github_local/DarwinFinches/Danny2024/D_CandidateGeneAnalysis/GeneLists/PCA

### High
R

library(ggvenn)

df_cra<-read.csv("cra.high_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("for.high_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("par.high_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

overlap_cra_for <- intersect(df_cra, df_for)
overlap_cra_par <- intersect(df_cra, df_par)
overlap_for_par <- intersect(df_for, df_par)
overlap_any_two <- union(overlap_cra_for, union(overlap_cra_par, overlap_for_par))
intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))

write(overlap_cra_for, file="overlap_cra_for.pca.high.all.csv")
write(overlap_cra_par, file="overlap_cra_par.pca.high.all.csv")
write(overlap_for_par, file="overlap_for_par.pca.high.all.csv")
write(overlap_any_two, file="overlap_any_two.pca.high.all.csv")
write(intersect_all, file="intersection.pca.high.all.csv")

pdf(file = "plots/intersection.pca.high.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()

### Highest

R

library(ggvenn)

df_cra<-read.csv("cra.highest_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("for.highest_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("par.highest_PC1_2.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

gene_sets <- list(
  cra= df_cra,
  "for" = df_for,
  par = df_par
)

overlap_cra_for <- intersect(df_cra, df_for)
overlap_cra_par <- intersect(df_cra, df_par)
overlap_for_par <- intersect(df_for, df_par)
overlap_any_two <- union(overlap_cra_for, union(overlap_cra_par, overlap_for_par))
intersect_all <- Reduce(intersect, list(df_cra, df_for, df_par))

write(overlap_cra_for, file="overlap_cra_for.pca.highest.all.csv")
write(overlap_cra_par, file="overlap_cra_par.pca.highest.all.csv")
write(overlap_for_par, file="overlap_for_par.pca.highest.all.csv")
write(overlap_any_two, file="overlap_any_two.pca.highest.all.csv")
write(intersect_all, file="intersection.pca.highest.all.csv")

pdf(file = "plots/intersection.pca.highest.all.pdf", width = 6, height = 6, useDingbats=FALSE)

ggvenn(gene_sets,
    columns = c("cra", "for", "par"),
  stroke_size = 0.5, set_name_size = 4
  )
dev.off()