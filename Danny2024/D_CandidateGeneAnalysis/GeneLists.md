# running interactively because it is a quick script
chmod +x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/genelist/
# FST
# Define species, environments, and window sizes
# Identify genes in windows
species=( "cra" "for" "par")
window_sizes=( 50000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
    -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${sp}pre_${sp}post/${sp}pre_${sp}post.fst_${win}.outlier.csv \
    -n ${sp} \
    -m fst \
    -w ${win}

    done
done

# Tajima diff

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

# chromosome 1

grep 'NC_044571.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr1.bed


GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr1.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr1.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr1.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr1.genenames.txt

# chromosome 9
grep 'NC_044579.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr9.bed

GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr9.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr9.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr9.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr9.genenames.txt

# chromosome 11
grep 'NC_044581.1' crapost.tajima_diff.50000kb.bed > crapost.tajima_diff.50000kb.chr11.bed

GFF=/xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/genomic.gff
bedtools intersect -a ${GFF} -b crapost.tajima_diff.50000kb.chr11.bed -wa > ../genes/crapost.tajima_diff.50000kb.chr11.genes.txt

grep 'ID\=gene' ../genes/crapost.tajima_diff.50000kb.chr11.genes.txt | awk '{OFS = "\t"} {split($9, arr, ";"); print(arr[1])}' | sed 's/ID\=gene\-//g' | sort -u > ../gene_names/crapost.tajima_diff.50000kb.chr11.genenames.txt



# raisd


species=( "cra_pre" "cra_post" "for_pre" "for_post" "par_pre" "par_post" )
window_sizes=( 50 )



# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    echo ${sp}
    ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/DarwinFinches/param_files/${sp}_params_raisd.sh \
    -i /xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${sp}/${sp}.raisd_${win}.outlier.csv \
    -n ${sp} \
    -m raisd \
    -w ${win}

    done
done

# RAiSD
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