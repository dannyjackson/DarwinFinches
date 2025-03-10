# running interactively because it is a quick script
chmod +x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh
# chmod -x ~/programs/DarwinFinches/Genomics-Main/general_scripts/bed_to_genelist.sh

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
