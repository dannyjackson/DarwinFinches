# C4_tajima.md
# zip some files up
sbatch --account=mcnew \
        --job-name=zipping \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.zipping.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=30:00:00 \
        --mem=50gb \
        /xdisk/mcnew/dannyjackson/cardinals/analyses/thetas/zipthemfilesup.sh

# run once per population to generate all files
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )

# Iterate over each combination
for sp in "${species[@]}"; do
sbatch --account=mcnew \
        --job-name=tajima_500000_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_500000_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=10:00:00 \
        --mem=50gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
        -w 500000 -s 500000
done




# Define species, environments, and window sizes
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )
window_sizes=( 50000 5000 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        step=$(expr $win / 4)

        sbatch --account=mcnew \
               --job-name=tajima_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.tajima_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=1:00:00 \
               --mem=50gb \
               ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
               -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
               -w $win -s $step
    done
done


# do snps
# Define species, environments, and window sizes
species=( "crapre" "crapost" "forpre" "forpost" "parpre" "parpost" )

# Iterate over each combination
for sp in "${species[@]}"; do

sbatch --account=mcnew \
        --job-name=tajima_${win}_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.tajima_${win}_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=1:00:00 \
        --mem=50gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/tajima/tajima.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_tajima.sh \
        -w 1 -s 1
done



# calculate change in Tajima
# 50,000 snps

# cra
(echo -e "chromo\tposition\tTajima"; awk 'BEGIN {OFS="\t"} NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' crapre/crapre.Tajima.50000.Ztransformed.csv crapost/crapost.Tajima.50000.Ztransformed.csv) | grep -v 'NW'  | grep -v 'Z' >> cradiff.txt


# for
(echo -e "chromo\tposition\tTajima"; awk 'NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' forpre/forpre.Tajima.50000.Ztransformed.csv forpost/forpost.Tajima.50000.Ztransformed.csv) >> fordiff.txt

# par
(echo -e "chromo\tposition\tTajima"; awk 'NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' parpre/parpre.Tajima.50000.Ztransformed.csv parpost/parpost.Tajima.50000.Ztransformed.csv) >> pardiff.txt

# plot these changes
# cra




# Repeat outlier and plotting analyses using depth and map filtered windows

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas

species=( "cra" "for" "par")
# first, add header to input file (double check if necessary before doing so)
for sp in "${species[@]}"; do

    PREFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}pre/50000/${sp}pre.theta.thetasWindow.pestPG.depthmapfiltered"
    
    POSTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}post/50000/${sp}post.theta.thetasWindow.pestPG.depthmapfiltered"

    (echo -e "chromo\tposition\tTajima"; awk 'BEGIN {OFS="\t"} NR==FNR && FNR>1 {data[$2,$3] = $9; next} FNR>1 && ($2,$3) in data {print $2, $3, $9 - data[$2,$3]}' ${PREFILE} ${POSTFILE}) >> ${sp}_tajimadiff.depthmapfiltered.txt

done


species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    source ~/programs/DarwinFinches/param_files/${sp}pre_params_tajima.sh

    FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${sp}_tajimadiff.depthmapfiltered.txt"

    cp ${FILE} "${FILE}.numchrom" 

    # Replace chromosome names if conversion file is provided
    if [ -n "$CHR_FILE" ]; then
        echo "Replacing chromosome names based on conversion file..."
        while IFS=',' read -r first second; do
            sed -i "s/$second/$first/g" "${FILE}.numchrom" 
        done < "$CHR_FILE"
    fi

    # z transform windowed data
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.filteredfiles.r" \
    "${OUTDIR}" "${CUTOFF}" "${FILE}.numchrom" 

    cp ${FILE} "${FILE}.numchrom" 

    Z_OUT="${FILE}.numchrom.Ztransformed.csv"

    # sed -i 's/\"//g' ${Z_OUT}

    # Run R script for plotting
    echo "Generating Manhattan plot from ${Z_OUT}..."
    Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.filteredfiles.tajimadiff.r" \
        "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${Z_OUT}" "Tajima" 

    echo "Script completed successfully!"

done