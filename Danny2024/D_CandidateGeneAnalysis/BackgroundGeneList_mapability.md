
cd /xdisk/mcnew/dannyjackson/cardinals/datafiles/snpable/
/xdisk/mcnew/dannyjackson/cardinals/datafiles/snpable/GCF_901933205_revised_mask.150.50.fa


#!/bin/bash

cd /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/

# Input file
FASTA_MASK="/xdisk/mcnew/dannyjackson/cardinals/datafiles/snpable/GCF_901933205_revised_mask.150.50.fa"
OUT_MASK="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_site_mapstats.txt"

awk '
  /^>/ {
    chrom = $1
    sub(/^>/, "", chrom)
    pos = 0
    next
  }
  {
    for (i = 1; i <= length($0); i++) {
      pos++
      base = substr($0, i, 1)
      print chrom, pos, (base == "3" ? 1 : 0)
    }
  }
' ${FASTA_MASK} > ${OUT_MASK}

# compute average mapability
grep 'NC_' "${OUT_MASK}"| grep -v "NC_044601" > "${OUT_MASK}.autosomes"


# Define file paths
DEPTH_FILE="${OUT_MASK}.autosomes"
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.bam"
AVGMAP_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.map.csv"


sbatch --account=mcnew \
--job-name=fstavgmap \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fstavgmap.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGMAP_FILE} -t 94
# 12606200

# repeat for thetas

DEPTH_FILE="${OUT_MASK}.autosomes"
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.bam"
AVGMAP_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.thetas.map.csv"



# run via generalized script on puma
sbatch --account=mcnew \
--job-name=thetasavgdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.thetasavgdepth.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGMAP_FILE} -t 94
# 12606350











# repeat with RAiSD
CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
OUT_MASK="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_site_mapstats.txt"
MAP_FILE="${OUT_MASK}.autosomes"

SCRIPT_PATH="/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh"
OUTDIR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd"
BAMSTATS_DIR="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats"
SLURM_OUT="slurm_output"

mkdir -p "${SLURM_OUT}"

for POP in cra_post cra_pre for_post for_pre par_post par_pre; do
# for POP in cra_pre; do

    WIN_FILE="${OUTDIR}/${POP}/50/RAiSD_${POP}.windows.bed"
    AVGMAP_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.map.csv"


    sbatch --account=mcnew \
        --job-name=fstavgmap_${POP} \
        --partition=standard \
        --mail-type=ALL \
        --output=${SLURM_OUT}/output.fstavgmap.${POP}.%j \
        --nodes=1 \
        --ntasks-per-node=94 \
        --time=48:00:00 \
        "${SCRIPT_PATH}" -d "${MAP_FILE}" -w "${WIN_FILE}" -a "${AVGMAP_FILE}" -t 94
done