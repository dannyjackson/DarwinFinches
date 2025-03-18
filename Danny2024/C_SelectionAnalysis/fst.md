# generate SAF files 
cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/safs

species=( "cra" "for" "par")

species=( "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=saf_${sp}_pre \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_pre.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}pre

    sbatch --account=mcnew \
            --job-name=saf_${sp}_post \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.saf_${sp}_post.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \
            -p /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh \
            -n ${sp}post
done

~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh 

/xdisk/mcnew/finches/dannyjackson/finches/referencelists/allsnps.sites_headless.mafs

Set the working directory to fst so that all relevant slurm output files appear here as well:
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst
```
# Run once per species to generate all files

# Define species, environments, and window sizes

# fst
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst

sp=( "cra" )
sp=( "for" )
sp=( "par" )

chmod +x ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh 

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
-p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
-w 500000 -s 500000

sp=( "par" )

sbatch --account=mcnew \
        --job-name=fst_500000_${sp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.fst_500000_${sp}.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=3:00:00 \
        --mem=100gb \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
        -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
        -w 500000 -s 500000
        
# Iterate over each combination
species=( "cra" "for" "par")

for sp in "${species[@]}"; do
    sbatch --account=mcnew \
            --job-name=fst_500000_${sp} \
            --partition=standard \
            --mail-type=ALL \
            --output=slurm_output/output.fst_500000_${sp}.%j \
            --nodes=1 \
            --ntasks-per-node=4 \
            --time=7:00:00 \
            --mem=100gb \
            ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
            -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
            -w 500000 -s 500000
done

~/programs/angsd/angsd -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doCounts 1 -doDepth 1 -setMinDepthInd 4 -minInd 20 -minQ 30 -minMapQ 30 -sites /xdisk/mcnew/dannyjackson/finches/angsty/analyses/sites_headless.mafs -bam /xdisk/mcnew/dannyjackson/finches/reference_lists/allsamplebams.txt -out /xdisk/mcnew/dannyjackson/finches/angsty/analyses/saf/allsnps -anc /xdisk/mcnew/dannyjackson/finch_wgs/fastqs/GCF_901933205.fa



# Iterate over several window sizes

# Define species, environments, and window sizes
species=( "cra" "for" "par")
window_sizes=( 50000 5000 )


# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do
        time_limit=${time_limits[$win]}
        step=$(expr $win / 4)
            ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
            -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh \
               -w $win -s $step 
    done
done



window_sizes=( 1 )

# Iterate over each combination
for win in "${window_sizes[@]}"; do
    for sp in "${species[@]}"; do

        sbatch --account=mcnew \
               --job-name=fst_${win}_${sp} \
               --partition=standard \
               --mail-type=ALL \
               --output=slurm_output/output.fst_${win}_${sp}.%j \
               --nodes=1 \
               --ntasks-per-node=4 \
               --time=1:00:00 \
               --mem=50gb \
                ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/fst/fst.sh \
                -p ~/programs/DarwinFinches/param_files/${sp}_params_fst.sh  \
               -w $win -s $win 
    done
done


