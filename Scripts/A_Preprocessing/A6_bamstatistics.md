# BAM Statistics & Depth Summaries — Annotated Guide
This computes statistics on bam files generated in previous preprocessing steps.
---

## Prerequisites
- Index BAMs (`.bai`) exist for any commands that need them (flagstat does not, depth does not strictly need it, but indexing is best practice).
- Input lists contain **sample IDs only** (no extensions), one per line.

---

## Flagstat & Depth (clipoverlap + sortedmarked) in parallel

> **What this does:** For each sample ID, computes `flagstat` and `depth` on both the *clipped* BAMs under `clipoverlap/` and the *sorted+marked* BAMs under `sortedmarkedbamfiles/`.

Create reference list:
```bash
cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/
ls *all* | awk 'BEGIN {FS = "."} {print $1}' > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt
```

Create batch script:
```bash
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=10:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

module load parallel
module load samtools 

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/

stating () { 
samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.clipped.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$@".depthstats

samtools flagstat /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam > /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".bamstats

samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".sorted.marked.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$@".depthstats

}

export -f stating 

parallel -j 12 stating :::: /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt  

```

Submit:
```bash
sbatch stating_dupmarkedbams.sh 
```


## Mean & SD of depth per sample (TSV output)

> **What this does:** Converts each `*.depthstats` (3 columns: `chrom pos depth`) into a single-line summary with `sample,mean,sd,n_sites`.

First, do it on sortedmarkedbamfiles.
Create sbatch:
```bash
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 

```
Submit
```bash
sbatch averagedepth_sortedbams.sh 
```
Then, do it on clipoverlap bams.
Create sbatch:
```bash
#!/bin/bash

#SBATCH --job-name=bamstats
#SBATCH --ntasks=12
#SBATCH --nodes=1             
#SBATCH --time=30:00:00   
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --mem-per-cpu=5gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dannyjackson@arizona.edu
#SBATCH --output=output.bamstats.%j

cd /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/


while read -r finch;
do 
  echo $finch
  echo $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

while read -r finch;
do 
  echo "all" $finch
  echo "all" $finch >> depthstats.txt

  # average and standard deviaiton
  awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' "$finch".all.depthstats >> depthstats.txt

done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_all_dupmarked.txt 

```

---

## Depth with mapping quality filter on realigned BAMs
```bash
while read -r finch;
do 
  samtools depth /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/sortedmarkedbamfiles/"$finch".sorted.marked.clipped.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/clipoverlap/"$finch".depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt 

# Compute per-sample sequencing depth for bases with mapping quality higher than 20


module load samtools
while read -r finch;
do 
  samtools depth -q 30 /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.bam -o /xdisk/mcnew/dannyjackson/finches/bias_testing/batchnaive/indelrealignment/"$finch".realigned.depthstats
done < /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_dupmarked.txt

# Create a proper list of realigned samples (names, not a count)

ls *bam | awk 'BEGIN {FS = "."} {print $1}' | wc -l > /xdisk/mcnew/dannyjackson/finches/reference_lists/samplenames_indelrealigned.txt
```

---
## ANGSD depthSample → sample-level average depth (Python)

> **What this does:** Interprets `stats.<CHROM>.depthSample` where columns represent coverage counts per depth (1..N), computes the **mean depth** per sample, and joins with a sample metadata table.


```python
python
import numpy as np
import pandas as pd 

# File where row i corresponds to sample i, and col j stores counts of sites with depth=j
# (ANGSD depthSample style)
df = pd.read_csv('stats.NC_044571.depthSample', sep='\t', header=None)
cols = df.shape[1]
multiplier = list(range(1, cols+1))

df_adj = df.mul(multiplier, axis = 1)

counts = df.sum(axis = 1)
totaldepth = df_adj.sum(axis = 1)
avgdepth = totaldepth / counts

# Join with metadata
df_names = pd.read_csv('/xdisk/mcnew/dannyjackson/finches/reference_lists/sample_species_treatment.txt', sep='\t')

df_final = pd.concat([df_names, avgdepth], axis=1)

df_final = df_final.rename({0: 'avgdepth'}, axis=1)

df_final.groupby(['species', 'treatment']).mean()

df_final.sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'CRA'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'FOR'].sort_values(by=['avgdepth']).round(decimals=2)

df_final[df_final["species"] == 'PAR'].sort_values(by=['avgdepth']).round(decimals=2)

print(df_final.groupby(['species', 'treatment']).size())

df_final.groupby(['species', 'treatment']).mean('MEAN_DEPTH')
```
---

## Quick one-liner: averages from 2‑column depth files
```bash
# compute average depth of individuals
for f in coverage_samtoolsDepth_*.txt; do
    avg=$(awk '{sum += $2} END {if (NR > 0) print sum / NR}' "$f")
    sample=$(basename "$f" .txt | sed 's/^coverage_samtoolsDepth_//')
    echo -e "${sample}\t${avg}"
done | sort -k2,2nr > indv_depthstats.txt
```