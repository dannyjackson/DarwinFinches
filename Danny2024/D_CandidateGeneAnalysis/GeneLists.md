# running interactively because it is a quick script
chmod +x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh

# FST
# Define species, environments, and window sizes
# Identify genes in windows
species=( "noca" "pyrr" )
window_sizes=( 50000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

    ~/programs/CardinalisGenomics/Genomics-Main/general_scripts/bed_to_genelist.sh \
    -p ~/programs/CardinalisGenomics/${sp}_params_fst.sh \
    -i /xdisk/mcnew/dannyjackson/cardinals/analyses/fst/${sp}urban_${sp}rural/${sp}urban_${sp}rural.fst_${win}.outlier.csv \
    -n ${sp} \
    -m fst \
    -w ${win}

    done
done
