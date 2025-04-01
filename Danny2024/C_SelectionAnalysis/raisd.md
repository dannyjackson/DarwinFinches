# RAiSD.md 

# merge vcfs
# path to lists of samples in each population

sed -i 's/PL/lamich_PL/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre_pops.txt
sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt
sed -i 's/PARV/lamich_PARV/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt
sed -i 's/SRR/SRR2917/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt


/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post_pops.txt 


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt
done

for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt
done


for IND in `cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre_pops.txt  `;
 do 
 bcftools concat /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}*gz >> /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf
 echo "/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/${IND}.allchrom.phased.vcf" >> /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt
done

for IND in `ls *vcf `
 do
 echo ${IND}
 bgzip ${IND}
 done

# FOR PRE
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

####### HERE #####


sbatch --account=mcnew \
--job-name=indexing \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.indexing.%j \
--nodes=1 \
--ntasks-per-node=12 \
--time=10:00:00 \
/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf2/indexing.sh 

#!/bin/bash

module load bcftools

# FOR POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done


# CRA PRE
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done


# CRA POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

# PAR PRE

sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

# PAR POST
sed -i 's/\.vcf/\.vcf\.gz/g' /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt

for IND in ` cat /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt `
 do
 echo ${IND}
 bcftools index ${IND}
 done

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/for_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/cra_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_pre.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf.gz

bcftools merge --file-list /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering/par_post.phasedvcfs.txt -Oz -o /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf.gz


#!/bin/bash

# Define arrays for species and time periods
species=("cra" "par" "for")
time_periods=("pre" "post")

# Base directory paths
file_list_base="/xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_filtering"
output_base="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3"

# Loop through species and time periods
for sp in "${species[@]}"; do
    for tp in "${time_periods[@]}"; do
        input_list="${file_list_base}/${sp}_${tp}.phasedvcfs.txt"
        output_file="${output_base}/${sp}_${tp}.phased.vcf.gz"

        sbatch --account=mcnew \
        --job-name=merge.${sp}.${tp} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/merge.${sp}.${tp}.%j \
        --nodes=1 \
        --ntasks-per-node=1 \
        --time=24:00:00 \
        merge.sh -i "$input_list" -o "$output_file"
        
    done
done

#!/bin/bash

while getopts ":i:o:" option; do
    case "${option}" in
        i) INPUT=${OPTARG} ;;
        o) OUTPUT=${OPTARG} ;;
        *) echo "Invalid option: -${OPTARG}" >&2; usage ;;
    esac
done

module load bcftools

bcftools merge --file-list "$INPUT" -Oz -o "$OUTPUT"

sbatch --account=mcnew \
    --job-name=raisd.test \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.raisd.test.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=10:00:00 \
    ~/programs/raisd_preprocessing.sh 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd

gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf.gz
gunzip /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf.gz

SPECIES=( "cra_pre" ) #  "cra_post" "for_pre" "for_post" "par_pre" "par_post" 
WINDOW_SIZES=( 50 ) # 10 20 50 100 500 1000 


# Iterate over each combination
for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    sbatch --account=mcnew \
        --job-name=raisd_${WIN}_${SP} \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.raisd_${WIN}_${SP}.%j \
        --nodes=1 \
        --time=1:00:00 \
        ~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/raisd/raisd.sh -p ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh -n ${SP} -w ${WIN}
    done 
done

~/programs/DarwinFinches/Genomics-Main/C_SelectionAnalysis/raisd/raisd.sh -p ~/programs/DarwinFinches/param_files/cra_pre_params_raisd.sh -n cra_pre -w 50 -r $REF

source /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_base.sh

SPECIES=( "cra_post" "for_pre" "for_post" "par_pre" "par_post" ) # "cra_pre" 
WINDOW_SIZES=( 50 ) # 10 20 50 100 500 1000 


for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do

    POP=${SP}
    source ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh
    
    cd ${OUTDIR}/analyses/raisd/${SP}/${WIN}

    ~/programs/RAiSD/raisd-master/RAiSD \
        -n "${SP}" \
        -I /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/"${SP}".phased.vcf \
        -f -O -R -P -C "${REF}" -w "${WIN}"
    
    mv *pdf plots
    mv *Report* reportfiles
    mv *Info* infofiles
    done
done


SP="cra_pre"
WIN="50"

SPECIES=( "for_pre") 

for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh

        echo -e "chrom\tposition\tstart\tend\tVAR\tSFS\tLD\tU" > ${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes


        while read -r SCAFFOLD;
        do 
        awk -v CHROM="$SCAFFOLD" '{print CHROM, $0}' ${OUTDIR}/analyses/raisd/${SP}/${WIN}/reportfiles/RAiSD_Report.${SP}.${SCAFFOLD} >> ${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes

        done < "${OUTDIR}/referencelists/SCAFFOLDS.txt"

        # rename scaffolds for plotting

        echo "Processing CHROM file: $CHR_FILE..."

        FILE="${OUTDIR}/analyses/raisd/${SP}/${WIN}/RAiSD_Report.${SP}.chromosomes"
        # Read CHROM line by line
        while IFS=',' read -r first second; do
            echo "Replacing occurrences of '$second' with '$first' in $FILE"
            sed -i.bak "s/$second/$first/g" "$FILE"
        done < "$CHR_FILE"

        rm -f "${FILE}.bak"

        sed -i 's/ /\t/g' "${FILE}"


        # z transform U metric


        echo 'z transforming u metric'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/ztransform_windows.r" "${OUTDIR}" "${CUTOFF}" "${FILE}" "${WIN}" "${SP}"
        echo 'finished z transforming' 

        WIN_OUT="${OUTDIR}/analyses/raisd/${SP}/${SP}.raisd.${WIN}.Ztransformed.csv"

        sed -i 's/\"//g' $WIN_OUT

        # plot scaffolds
        echo 'visualizing windows'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "raisd" "${SP}"
        echo 'finished windowed plot' 
    done
done

for WIN in "${WINDOW_SIZES[@]}"; do
    for SP in "${SPECIES[@]}"; do
        source  ~/programs/DarwinFinches/param_files/${SP}_params_raisd.sh
        # plot scaffolds
        echo 'visualizing windows'
        Rscript "${SCRIPTDIR}/Genomics-Main/general_scripts/manhattanplot.r" "${OUTDIR}" "${COLOR1}" "${COLOR2}" "${CUTOFF}" "${WIN_OUT}" "${WIN}" "raisd" "${SP}"
        echo 'finished windowed plot' 
    done
done