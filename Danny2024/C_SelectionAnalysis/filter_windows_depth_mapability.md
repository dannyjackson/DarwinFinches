# Filter windows by depth and mapability .md

After running the windowed analyses but BEFORE computing outliers or plotting the output, I filtered out windows with low or high average statistic across the window. For instance, windows with low average depth may drive unusual patterns and erroneously appear as an outlier. I filtered for average depth and average mapability (see scripts related to snpable), keeping only windows within 2 standard deviations of the mean for both statistics across all analyses.

### Compute average statistic per site across the genome. Include every site from BAM files, not just SNPs

#### depth
```
#!/bin/bash

awk 'NR==1 {next} { 
    sum = 0; 
    for (i = 3; i <= NF; i++) sum += $i; 
    avg = sum / (NF - 2); 
    print $1, $2, avg;
}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_depthstats.txt  > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt 


sbatch --account=mcnew \
--job-name=avgdepthstats \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.avgdepthstats.%j \
--nodes=1 \
--ntasks-per-node=1 \
--time=5:00:00 \
/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/avgdepthstats.sh

```

#### mapability

```
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
```

# Make a bed file from the windowed analysis. Then run statavg_over_bedwindows.sh to compute average statistic per window


#### FST
```
awk 'BEGIN {OFS="\t"} {print $2, ($3-25000), $3+2500}' /xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/crapre_crapost/50000/crapre_crapost.50000.fst.chrom.autosomes > /xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.bam


# Define file paths
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt"

# define bam file
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.bam"

# define avg depth output file
AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/50kbwin.fst.depth.csv"


# run via generalized script on puma
sbatch --account=mcnew \
--job-name=fstavgdepthfastparallel3 \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fstavgdepthfastparallel3.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGDEPTH_FILE} -t 94
```

#### Thetas (Tajima's D)
```
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/chrom_avg_depthstats/chrom_avg_depthstats.txt"
WIN_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/windowed_bamfiles/50kbwin.thetas.bam"
AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.thetas.depth.csv"

awk 'BEGIN { OFS="\t" } NR>1 {print $2, ($3-25000), ($3+25000)}' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/50000/crapost.theta.thetasWindow.pestPG > "${WIN_FILE}.numchrom"


# Build the replacement logic into an AWK program
awk -v chr_file="$CHR_FILE" '
BEGIN {
    FS = OFS = "\t"
    while ((getline < chr_file) > 0) {
        split($0, a, ",")
        map[a[1]] = a[2]
    }
}
{
    if ($1 in map) $1 = map[$1]
    print
}
' "${WIN_FILE}.numchrom" | grep "NC_" | grep -v ${SEXCHR} > "${WIN_FILE}" 


sbatch --account=mcnew \
--job-name=thetasavgdepth \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.fstavgdepth.%j \
--nodes=1 \
--ntasks-per-node=94 \
--time=24:00:00 \
/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh -d ${DEPTH_FILE} -w ${WIN_FILE} -a ${AVGDEPTH_FILE} -t 94
```

#### RAiSD
```
CHR_FILE="/xdisk/mcnew/dannyjackson/cardinals/referencelists/GCF_901933205_chromconversion.txt"
DEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/chrom_avg_depthstats.txt"
SCRIPT_PATH="/home/u15/dannyjackson/programs/DarwinFinches/Genomics-Main/general_scripts/statavg_over_bedwindows.sh"
OUTDIR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd"
BAMSTATS_DIR="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats"
SLURM_OUT="slurm_output"

mkdir -p "${SLURM_OUT}"

# for POP in cra_post cra_pre for_post for_pre par_post par_pre; do

for POP in cra_pre; do
    RAISD_FILE="${OUTDIR}/${POP}/50/RAiSD_Report.${POP}.chromosomes"
    WIN_FILE="${OUTDIR}/${POP}/50/RAiSD_${POP}.windows.bed"
    AVGDEPTH_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.depth.csv"

    awk '{print $1, $3, $4}' "$RAISD_FILE" > "${WIN_FILE}.numchrom"

    awk -v chr_file="$CHR_FILE" '
    BEGIN {
        FS = OFS = " "
        while ((getline < chr_file) > 0) {
            split($0, a, ",")
            map[a[1]] = a[2]
        }
    }
    {
        if ($1 in map) $1 = map[$1]
        print
    }
    ' "${WIN_FILE}.numchrom" | grep 'NC_' | grep -v "NC_044601" > "${WIN_FILE}"

    sbatch --account=mcnew \
        --job-name=fstavgdepth_${POP} \
        --partition=standard \
        --mail-type=ALL \
        --output=${SLURM_OUT}/output.fstavgdepth.${POP}.%j \
        --nodes=1 \
        --ntasks-per-node=94 \
        --time=24:00:00 \
        "${SCRIPT_PATH}" -d "${DEPTH_FILE}" -w "${WIN_FILE}" -a "${AVGDEPTH_FILE}" -t 94
done
```

# Run filter_statavg_output.r (short script, easily run interactively)
### Depth
#### FST 
```
source ~/programs/DarwinFinches/param_files/params_base.sh 

AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.fst.depth.III.csv"

OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.fst.depth.filtered.bam"

grep "${CHR_LEAD}" ${AVGDEPTH_FILE} | grep -v "${SEXCHR}" > "${AVGDEPTH_FILE}.autosomes"

Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGDEPTH_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

rm "${OUTPUT_FILE}.unsorted"

# Filter FST files
BAMFILE="${OUTPUT_FILE}"
awk '{print $1, ($2+25000)}' "${BAMFILE}" > "${BAMFILE}.midpos"


for POP in cra for par; do
    FSTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${POP}pre_${POP}post/50000/${POP}pre_${POP}post.50000.fst"
    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${POP}pre_${POP}post/50000/${POP}pre_${POP}post.50000.fst.depthfiltered"
    grep "${CHR_LEAD}" ${FSTFILE} | grep -v "${SEXCHR}" > "${FSTFILE}.autosomes"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in FSTFILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $2
        midpoint = $3
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "$BAMFILE" "$FSTFILE" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

    
done

```
#### Tajima
```
source ~/programs/DarwinFinches/param_files/params_base.sh 


AVGDEPTH_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.thetas.depth.csv"

OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.thetas.depth.filtered.bam"
BAMFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.thetas.depth.filtered.bam.midpos"

grep "${CHR_LEAD}" ${AVGDEPTH_FILE} | grep -v "${SEXCHR}" > "${AVGDEPTH_FILE}.autosomes"

Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGDEPTH_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

rm "${OUTPUT_FILE}.unsorted"

awk '{print $1, ($2+25000)}' "${OUTPUT_FILE}" > "${BAMFILE}"

# Filter thetas file

for POP in crapre crapost forpre forpost parpre parpost; do
  
    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${POP}/50000/${POP}.theta.thetasWindow.pestPG.depthfiltered"

    THETAFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${POP}/50000/${POP}.theta.thetasWindow.pestPG"


    awk -v chr_file="$CHR_FILE" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline < chr_file) > 0) {
            split($0, a, ",")
            map[a[1]] = a[2]
        }
    }
    {
        if ($2 in map) $2 = map[$2]
        print
    }
    ' "${THETAFILE}" | grep 'NC_' > "${THETAFILE}.txtchrom" 

    grep "${CHR_LEAD}" ${THETAFILE}.txtchrom | grep -v "${SEXCHR}" > "${THETAFILE}.autosomes"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in THETAFILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $2
        midpoint = $3
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "$BAMFILE" "$THETAFILE.autosomes" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

done
```


#### RAiSD
```
# cra_post for_post for_pre par_post par_pre

source ~/programs/DarwinFinches/param_files/params_base.sh 

for POP in cra_post cra_pre for_post for_pre par_post par_pre; do

for POP in cra_pre; do

    BAMSTATS_DIR="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed"
    OUTDIR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd"
    AVGDEPTH_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.depth.csv"
    RAISD_FILE="${OUTDIR}/${POP}/50/RAiSD_Report.${POP}.chromosomes"
    OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/depthstats/avgdepth_windowed/50kbwin.raisd.${POP}.depth.filtered.bam"

    
    grep "${CHR_LEAD}" ${AVGDEPTH_FILE} | grep -v "${SEXCHR}" > "${AVGDEPTH_FILE}.autosomes"

    Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGDEPTH_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

    sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    rm "${OUTPUT_FILE}.unsorted"

    # Filter RAiSD files
    BAMFILE="${OUTPUT_FILE}"
    awk '{print $1, (($2+$3)/2)}' "${BAMFILE}" > "${BAMFILE}.midpos"


    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${POP}/50/${POP}.50.raisd.depthfiltered"


    awk -v chr_file="$CHR_FILE" '
    BEGIN {
        FS = OFS = "\t"
        while ((getline < chr_file) > 0) {
            split($0, a, ",")
            map[a[1]] = a[2]
        }
    }
    {
        if ($1 in map) $1 = map[$1]
        print
    }
    ' "${RAISD_FILE}" | grep 'NC_' > "${RAISD_FILE}.txtchrom" 

    grep "${CHR_LEAD}" ${RAISD_FILE}.txtchrom | grep -v "${SEXCHR}" > "${RAISD_FILE}.autosomes"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in RAISD_FILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $1
        midpoint = $2
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "${BAMFILE}.midpos" "${RAISD_FILE}.autosomes" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

        
done

```
### Mapability
#### FST 
```
source ~/programs/DarwinFinches/param_files/params_base.sh 

AVGMAP_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.fst.map.csv"

OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.fst.map.filtered.bam"

grep "${CHR_LEAD}" ${AVGMAP_FILE} | grep -v "${SEXCHR}" > "${AVGMAP_FILE}.autosomes"

Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGMAP_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

rm "${OUTPUT_FILE}.unsorted"

# Filter FST files
BAMFILE="${OUTPUT_FILE}"
awk '{print $1, ($2+25000)}' "${BAMFILE}" > "${BAMFILE}.midpos"

# PAR stats
# depth filtered alone: 75010
# mapfiltered alone: 74110
# mapdepth: 73887

for POP in cra for par; do
    FSTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${POP}pre_${POP}post/50000/${POP}pre_${POP}post.50000.fst.depthfiltered"
    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/fst/${POP}pre_${POP}post/50000/${POP}pre_${POP}post.50000.fst.depthmapfiltered"
    grep "${CHR_LEAD}" ${FSTFILE} | grep -v "${SEXCHR}" > "${FSTFILE}.autosomes"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in FSTFILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $2
        midpoint = $3
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "$BAMFILE" "$FSTFILE" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

    
done

```
#### Tajima
```
source ~/programs/DarwinFinches/param_files/params_base.sh 


AVGMAP_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.thetas.map.csv"

OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.thetas.map.filtered.bam"
BAMFILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.thetas.map.filtered.bam.midpos"

grep "${CHR_LEAD}" ${AVGMAP_FILE} | grep -v "${SEXCHR}" > "${AVGMAP_FILE}.autosomes"

Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGMAP_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

rm "${OUTPUT_FILE}.unsorted"

awk '{print $1, ($2+25000)}' "${OUTPUT_FILE}" > "${BAMFILE}"


# Filter thetas file

for POP in crapre crapost forpre forpost parpre parpost; do
    THETAFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${POP}/50000/${POP}.theta.thetasWindow.pestPG.depthfiltered"

    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/${POP}/50000/${POP}.theta.thetasWindow.pestPG.depthmapfiltered"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in THETAFILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $2
        midpoint = $3
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "$BAMFILE" "$THETAFILE" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

done
```


#### RAiSD
```

source ~/programs/DarwinFinches/param_files/params_base.sh 

for POP in cra_post cra_pre for_post for_pre par_post par_pre; do

    BAMSTATS_DIR="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed"
    OUTDIR="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd"
    AVGMAP_FILE="${BAMSTATS_DIR}/50kbwin.raisd.${POP}.map.csv"
    RAISD_FILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${POP}/50/${POP}.50.raisd.depthfiltered"
    OUTPUT_FILE="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/bamstats/mapping_stats/avgmap_windowed/50kbwin.raisd.${POP}.map.filtered.bam"

    
    grep "${CHR_LEAD}" ${AVGMAP_FILE} | grep -v "${SEXCHR}" > "${AVGMAP_FILE}.autosomes"

    Rscript ~/programs/DarwinFinches/Genomics-Main/general_scripts/filter_statavg_output.r "${AVGMAP_FILE}.autosomes" "${OUTPUT_FILE}.unsorted"

    sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    rm "${OUTPUT_FILE}.unsorted"

    # Filter RAiSD files
    BAMFILE="${OUTPUT_FILE}"
    awk '{print $1, (($2+$3)/2)}' "${BAMFILE}" > "${BAMFILE}.midpos"


    OUTFILE="/xdisk/mcnew/finches/dannyjackson/finches/analyses/raisd/${POP}/50/${POP}.50.raisd.depthmapfiltered"

    awk '
    # Load BAMFILE: store all (chrom, midpoint) pairs
    FNR==NR {
        bam[$1, $2] = 1
        next
    }

    # For each line in RAISD_FILE, keep only if the (chrom, midpoint) pair is in BAMFILE
    {
        chrom = $1
        midpoint = $2
        if ((chrom, midpoint) in bam) {
            print
        }
    }
    ' "${BAMFILE}.midpos" "${RAISD_FILE}" > "$OUTFILE"

    # sort -k1,1 -k2,2n "${OUTPUT_FILE}.unsorted" > "${OUTPUT_FILE}"

    # rm "${OUTPUT_FILE}.unsorted"

        
done