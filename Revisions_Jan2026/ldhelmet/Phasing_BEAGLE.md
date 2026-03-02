# Phase With Beagle
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/ldhelmet/phased_genomes

## Create VCFs for each species
### CRA
sbatch --account=mcnew \
    --job-name=vcf_cra \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.vcf_cra.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_AllPop_VCF.sh \
    -p ~/programs/DarwinFinches/params_preprocessing.sh \
    -b /xdisk/mcnew/finches/dannyjackson/finches/referencelists/cra_all_bams.txt \
    -v /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/cra_all.vcf.gz 

### FOR
sbatch --account=mcnew \
    --job-name=vcf_for \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.vcf_for.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_AllPop_VCF.sh \
    -p ~/programs/DarwinFinches/params_preprocessing.sh \
    -b /xdisk/mcnew/finches/dannyjackson/finches/referencelists/for_all_bams.txt \
    -v /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/for_all.vcf.gz 

### PAR
sbatch --account=mcnew \
    --job-name=vcf_par \
    --partition=standard \
    --mail-type=ALL \
    --output=slurm_output/output.vcf_par.%j \
    --nodes=1 \
    --ntasks-per-node=8 \
    --time=48:00:00 \
    ~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A2.5_AllPop_VCF.sh \
    -p ~/programs/DarwinFinches/params_preprocessing.sh \
    -b /xdisk/mcnew/finches/dannyjackson/finches/referencelists/par_all_bams.txt \
    -v /xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/par_all.vcf.gz 

## Phase each VCF with Beagle5

### Create mamba environment
```
module load micromamba

micromamba create -n beagle \
  -c conda-forge -c bioconda \
  beagle=5.4 \
  bcftools \
  htslib \
  openjdk
```
## Run beagle phasing
```
micromamba activate beagle

# CRA
beagle \
    gt=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/cra_all.vcf.gz \
    out=cra_phased

# FOR
beagle -Xmx8g \
    gt=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/for_all.vcf.gz \
    out=for_phased

# PAR
beagle \
    gt=/xdisk/mcnew/finches/dannyjackson/finches/datafiles/vcf_for_beagle/par_all.vcf.gz \
    out=par_phased


# Basic stats (replacement for whatshap stats)
tabix cra_phased.vcf.gz
bcftools stats cra_phased.vcf.gz > crastats.txt

# Quick sanity check: do we now have phased genotypes?

bcftools view cra_phased.vcf.gz \
| bcftools query  -f '[%GT\n]' 2>/dev/null \
| awk '{tot++; if($0~"\\|") ph++; if($0~"/") un++} END{print "total",tot,"phased(|)",ph+0,"unphased(/)",un+0}'
