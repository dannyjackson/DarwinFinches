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
tail -n +12 par.fst.50000kb.panther.txt > par.fst.50000kb.panther.clean.txt

sed -i '' 's/par.fst.50000kb.genenames.txt //g' par.fst.50000kb.panther.clean.txt

sed -i '' 's/(raw P-value)/P-value/g' par.fst.50000kb.panther.clean.txt

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r par.fst.50000kb.panther.clean

head -1 par.fst.50000kb.panther.clean.txt > par.fst.50000kb.panther.clean.p01.txt

awk 'BEGIN {FS = "\t"} $7 < 0.01' par.fst.50000kb.panther.clean.txt >> par.fst.50000kb.panther.clean.p01.txt

# grep -o 'GO:[0-9]\+' cra.fst.50000kb.panther.clean.txt | head > go_terms.txt


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

  # bottom
  tail -n +12 ${sp}.tajima_diff.50000kb.bottom.panther.txt > ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  sed -i '' "s/${sp}.tajima_diff.50000kb.bottom.genenames.txt //g" ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  sed -i '' "s/(raw P-value)/P-value/g" ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt

  Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/panther_fdr.r ${sp}.tajima_diff.50000kb.bottom.panther.clean

  head -1 ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt > ${sp}.tajima_diff.50000kb.bottom.panther.clean.p01.txt

  awk 'BEGIN {FS = "\t"} $7 < 0.01' ${sp}.tajima_diff.50000kb.bottom.panther.clean.txt >> ${sp}.tajima_diff.50000kb.bottom.panther.clean.p01.txt

  # top

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



# Iterate over each combination
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

# make background gene lists

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


# filter background gene lists locally
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


# for unique genes in cra post that were not in cra pre
# 1. Genes in df_cra_post that are not in df_cra_pre
df_cra_genes<-read.csv("cra_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_cra_background<-read.csv("background_genes.cra_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_cra_post <- intersect(df_cra_genes, df_cra_background)

write(unique_cra_post, file="cra.raisd.txt")

# 2. Genes in df_for_post that are not in df_for_pre
df_for_genes<-read.csv("for_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for_background<-read.csv("background_genes.for_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_for_post <- intersect(df_for_genes, df_for_background)

write(unique_for_post, file="for.raisd.txt")


# 3. Genes in df_par_post that are not in df_par_pre
df_par_genes<-read.csv("par_post.raisd.50kb.unique.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par_background<-read.csv("background_genes.par_pre.raisd.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

unique_par_post <- intersect(df_par_genes, df_par_background)

write(unique_par_post, file="par.raisd.txt")

# evaluate intersection of all raisd genes

cra_for_raisd_intersect <- intersect(unique_cra_post, unique_for_post)
cra_par_raisd_intersect <- intersect(unique_cra_post, unique_par_post)
par_for_raisd_intersect <- intersect(unique_par_post, unique_for_post)

# intersection of all background genes
two_raisd_intersect_background <- intersect(cra_intersect, for_intersect)
all_raisd_intersect_background <- intersect(two_raisd_intersect_background, par_intersect)
write(all_raisd_intersect_background, file="background.all.raisd.txt")

# intersection of all candidate genes
all_raisd_intersect <- intersect(cra_for_raisd_intersect, unique_par_post)
filtered_par_intersect <- intersect(all_raisd_intersect, all_raisd_intersect_background)
write(all_raisd_intersect, file="all.raisd.txt")

# visualize overlapping genes

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

# Make plots of FDR passing genes
# RAISD
##### Only CRA had any significant GO terms
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.raisd.panther.clean.fdr.txt > cra.raisd.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.raisd.panther.clean.fdr.GO

# FST
# cra.fst.50000kb.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.fst.50000kb.panther.clean.fdr.txt > cra.fst.50000kb.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fst.50000kb.panther.clean.fdr.GO

# for.fst.50000kb.panther.clean.fdr.txt -- only has one GO term
grep -o 'GO:[0-9]\{7\}' pantheroutput/for.fst.50000kb.panther.clean.fdr.txt > for.fst.50000kb.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.fst.50000kb.panther.clean.fdr.GO


# Tajima
# cra.tajima_diff.50000kb.bottom.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.tajima_diff.50000kb.bottom.panther.clean.fdr.txt > cra.tajima_diff.50000kb.bottom.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.tajima_diff.50000kb.bottom.panther.clean.fdr.GO

# cra.tajima_diff.50000kb.top.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/cra.tajima_diff.50000kb.top.panther.clean.fdr.txt > cra.tajima_diff.50000kb.top.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.tajima_diff.50000kb.top.panther.clean.fdr.GO

# for.tajima_diff.50000kb.top.panther.clean.fdr.txt
grep -o 'GO:[0-9]\{7\}' pantheroutput/for.tajima_diff.50000kb.top.panther.clean.fdr.txt > for.tajima_diff.50000kb.top.panther.clean.fdr.GO.txt
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.tajima_diff.50000kb.top.panther.clean.fdr.GO





# generate all gene lists 
cat cra* | sort -u > cra.fstraisd.genenames.txt
cat for* | sort -u > for.fstraisd.genenames.txt
cat par* | sort -u > par.fstraisd.genenames.txt

# prepare GO term lists for visualization
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


# visualize
# Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fst.50000kb.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra_post.raisd.50kb.unique.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r cra.fstraisd.GO

# Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for.fst.50000kb.genenames.fdr
# Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r for_post.raisd.50kb.unique.genenames.fdr

Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r par.fst.50000kb.genenames.fdr
Rscript /Users/danjack/Documents/Github_local/Genomics-Main/D_GeneVisualization/D2_visualization.local.r par_post.raisd.50kb.unique.genenames.fdr
















# Analyze overlap between methods
library(ggvenn)

df_cra<-read.csv("GeneLists/cra.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_for<-read.csv("GeneLists/for.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]
df_par<-read.csv("GeneLists/par.fstraisd.genenames.txt", header=FALSE, stringsAsFactors = FALSE)[[1]]

# compare with owl and great tit
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



# test for significance in overlap
# cra vs for
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

# cra vs par
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

# for vs par
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