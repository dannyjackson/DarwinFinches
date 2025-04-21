# allele frequency change md 


# use SAF files generated during FST analyses
# /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/deltaAF

# cra pre
~/programs/angsd/angsd -b /xdisk/mcnew/finches/dannyjackson/finches/referencelists/craprebams.txt \
  -ref /xdisk/mcnew/dannyjackson/cardinals/datafiles/referencegenome/ncbi_dataset/data/GCF_901933205.1/GCF_901933205.1_STF_HiC_genomic.fna \
  -out cra_pre_postmafs \
  -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 \
  -minMapQ 30 -minQ 20 \
  -GL 1 \
  -doMaf 1 -doMajorMinor 1 -doPost 1 \
  -pest /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs/crapre.sfs \
  -minInd 3


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