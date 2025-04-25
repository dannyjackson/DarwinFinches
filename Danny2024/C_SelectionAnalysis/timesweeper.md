# script for timesweeper

### must be run on elgato -- breaks with slim > 4
### NOTE: the timesweeper code overwrites the output csv file so you have to edit it to read as follows:
#### file to edit: /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/lib/python3.8/site-packages/timesweeper/find_sweeps_vcf.py
#### predictions.to_csv(outfile, mode="a", header=not os.path.exists(outfile), index=False, sep="\t")
### Compute statistics that are used as parameters in timesweeper
```
# compute average missingness from vcf
module load vcftools
vcftools --vcf /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.vcf --missing-indv --out /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic

awk 'NR > 1 {sum += $5; count++} END {if (count > 0) print sum/count}' /xdisk/mcnew/finches/dannyjackson/finches/datafiles/genotype_calls/finches_snps_multiallelic.imiss
### 0.0206774

# computed Ne in a different script, see B_PopulationStructureAnalysis/EffectivePopulationSize.md
```
### Load in environment (hashtag is command used to intially create it)
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper

interactive -a mcnew -t 4:00:00 

module load python/3.8/3.8.12
# python3 -m venv --system-site-packages timesweeper_env
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 
```

### cra
#### merge vcf for input
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"
```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out}" | bgzip > "${vcf_in}"
```
#### Run timesweeper
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

timesweeper sim_custom -y cra.yaml
timesweeper condense -o cra.pkl -m 0.02 -y  cra.yaml --hft
timesweeper train -i cra.pkl -y cra.yaml --hft
# Mean absolute error for Sel Coeff predictions: 0.12233152240514755

#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra

vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/cra_all.phased.sorted.cleaned.vcf.gz"
timesweeper detect -i "${vcf_in}" -o /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/cra/cra_detect -y cra.yaml --hft

sbatch --account=mcnew \
    --job-name=detect_cra \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_cra.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_cra.sh
    # Submitted batch job 3938157
```


### for
#### merge vcf for input
```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/data`files/vcf3/for_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"
```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out}" | bgzip > "${vcf_in}"
```
#### Run timesweeper
```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

timesweeper sim_custom -y for.yaml
timesweeper condense -o for.pkl -m 0.02 -y  for.yaml --hft
timesweeper train -i for.pkl -y for.yaml --hft
# Mean absolute error for Sel Coeff predictions: 0.09318568557500839

#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for

vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/for_all.phased.sorted.cleaned.vcf.gz"
timesweeper detect -i "${vcf_in}" -o /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/for/for_detect -y for.yaml --hft

sbatch --account=mcnew \
    --job-name=detect_for \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_for.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_for.sh
    # Submitted batch job 3938106
```
### par
#### merge vcf for input
```
module load bcftools htslib

vcf_pre="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.vcf"
vcf_pre_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_pre.phased.sorted.vcf.gz"
bgzip ${vcf_pre}
bcftools sort "${vcf_pre}.gz" -Oz -o  "${vcf_pre_sorted}" 
bcftools index -t "${vcf_pre_sorted}" 

vcf_post="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.vcf"
vcf_post_sorted="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_post.phased.sorted.vcf.gz"

bgzip ${vcf_post}
bcftools sort "${vcf_post}.gz" -Oz -o  "${vcf_post_sorted}" 
bcftools index -t "${vcf_post_sorted}" 

vcf_out="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.vcf.gz"

bcftools merge "${vcf_pre_sorted}" "${vcf_post_sorted}" -Oz -o "${vcf_out}"
```
#### clean vcf
```
vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.cleaned.vcf.gz"

zgrep -v 'ID=VDB,.*Version=' "${vcf_out}" | bgzip > "${vcf_in}"
```
#### Run timesweeper
```

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par

timesweeper sim_custom -y par.yaml
timesweeper condense -o par.pkl -m 0.02 -y  par.yaml --hft
timesweeper train -i par.pkl -y par.yaml --hft
# Mean absolute error for Sel Coeff predictions: 0.02001003921031952

#!/bin/bash

module load python/3.8/3.8.12
source /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/timesweeper_env/bin/activate
module load slim/3.7.1 samtools bcftools 

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par

vcf_in="/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf3/par_all.phased.sorted.cleaned.vcf.gz"
timesweeper detect -i "${vcf_in}" -o /xdisk/mcnew/finches/dannyjackson/finches/analyses/timesweeper/analyses/par/par_detect -y par.yaml --hft

sbatch --account=mcnew \
    --job-name=detect_par \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.detect_par.%j \
    --nodes=1 \
    --ntasks-per-node=1 \
    --time=7:00:00 \
    detect_par.sh
    # Submitted batch job 3938105
```