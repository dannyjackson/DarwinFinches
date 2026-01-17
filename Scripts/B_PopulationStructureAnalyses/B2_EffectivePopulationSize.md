cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop

# This computes thetas:
sbatch --array=0-5 B2_1_run_thetas_array.sh
sbatch --array=0-5 B2_2_run_thetas_array.sh
sbatch --array=0 B2_2_run_thetas_array.sh # 19071514

# Then I'll have to run something like this:

module load micromamba 
micromamba activate r_elgato

#!/bin/bash

POPS=(crapre crapost forpre forpost parpre parpost)

for POP in "${POPS[@]}"; do
  echo "Processing ${POP}"

  IN="thetas/${POP}.thetas.idx.pestPG"
  OUT="thetas/${POP}.thetas.idx.pestPG.autosomes"

  # write header
  head -n 1 "$IN" > "$OUT"

  # append autosomes (keep NC_, drop Z = NC_044601)
  grep 'NC_' "$IN" | grep -v 'NC_044601' >> "$OUT"
done

pops <- c("crapre", "crapost", "forpre", "forpost", "parpre", "parpost")

pops <- c("crapre", "forpre", "forpost", "parpre", "parpost")

mu <- 2.04e-9

results <- data.frame(pop=character(), theta_per_site=numeric(), Ne=numeric(), stringsAsFactors=FALSE)

for (pop in pops) {
  file <- paste0("thetas/", pop, ".thetas.idx.pestPG.autosomes")

  df <- read.table(
    file,
    header = TRUE,
    sep = "",            # any whitespace
    comment.char = "",   # IMPORTANT: keep the '#...' header line
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  df$tW     <- as.numeric(df$tW)
  df$nSites <- as.numeric(df$nSites)

  theta_per_site <- sum(df$tW, na.rm = TRUE) / sum(df$nSites, na.rm = TRUE)
  Ne <- theta_per_site / (4 * mu)

  results <- rbind(results, data.frame(pop = pop, theta_per_site = theta_per_site, Ne = Ne))
}












θ = 2N e(f)μ, where μ is the substitution rate per generation with estimations of thet

# From Lamichhaney et al. 2015:
# We used the following previously reported estimated mutation rates for nuclear and mtDNA: nuclear DNA, 2.04 × 10−9 per site per year estimated from the synonymous mutation rate on the Darwin’s finches’ lineage since the split from zebra finch45; mtDNA, a fossil-calibrated divergence rate of 2.1% per million years for bird cytochrome b sequences46.

# first, compute thetas using all sites:
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop
~/programs/DarwinFinches/Genomics-Main/A_Preprocessing/A1.3_siteallelefrequency.sh \


#!/bin/bash
source /home/u15/dannyjackson/programs/DarwinFinches/param_files/params_preprocessing.sh

POP=crapre
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches

cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop

~/programs/angsd/angsd \
  -bam ${OUTDIR}/referencelists/${POP}bams.txt \
  -rf ${OUTDIR}/referencelists/autosomes.txt \
  -out /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop/safs/${POP} \
  -dosaf 1 \
  -GL 1 \
  -doGlf 2 \
  -doMajorMinor 1 \
  -doCounts 1 \
  -doDepth 1 \
  -setMinDepthInd ${MINDEPTHIND} \
  -minInd ${MININD} \
  -minQ ${MINQ} \
  -minMapQ ${MINMAPQ} \
  -anc ${REF} \
  -nThreads 8

sbatch --account=mcnew \
        --job-name=saf_crapre \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.saf_crapre.%j \
        --nodes=1 \
        --ntasks-per-node=8 \
        --time=7:00:00 \
        --mem=100gb \
        crapre_SAFs.sh


#!/bin/bash
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches
POP=crapre
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop
~/programs/angsd/misc/realSFS -P 24 "/xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop/safs/${POP}.saf.idx" > /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop/sfs_files/crapre.sfs

sbatch --account=mcnew \
        --job-name=crapre_sfs \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.crapre_sfs.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=2:00:00 \
        --mem=50gb \
        crapre_sfs.sh

# Compute per-site thetas
#!/bin/bash
OUTDIR=/xdisk/mcnew/finches/dannyjackson/finches
POP=crapre
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop
~/programs/angsd/misc/realSFS saf2theta "${OUTDIR}/datafiles/safs/${POP}.saf.idx" \
    -outname "/xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop/thetas/${POP}" \
    -sfs crapre.sfs

sbatch --account=mcnew \
        --job-name=crapre_thetas \
        --partition=standard \
        --mail-type=ALL \
        --output=slurm_output/output.crapre_thetas.%j \
        --nodes=1 \
        --ntasks-per-node=4 \
        --time=2:00:00 \
        --mem=50gb \
        crapre_thetas.sh

# Compute average Watterson's theta for each pop
# cra pre


realSFS saf.idx -P 8 > crapre.sfs
realSFS saf.idx -P 8 -fold 1 > crapre.folded.sfs
realSFS saf.idx -P 8 -theta crapre.theta

~/programs/angsd_elgato/angsd/misc/thetaStat do_stat thetas/crapre.thetas.idx


head -n 1 thetas/crapre.thetas.idx.pestPG > thetas/crapre.thetas.idx.pestPG.autosomes
grep 'NC_' thetas/crapre.thetas.idx.pestPG | grep -v 'NC_044601' >> thetas/crapre.thetas.idx.pestPG.autosomes


df <- read.csv("thetas/crapre.thetas.idx.pestPG.autosomes", sep='\t')
theta_per_site <- sum(df$tW) / sum(df$nSites)
Ne <- theta_per_site / (4*2.04e-9)

155855.2



# cra post
head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819
tW_per_site / (4 * 2.04e-9)
## 90,730.43
## mean: 95354.65

# for pre
head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819

tW_per_site / (4 * 2.04e-9)
## 206,862.1

# for post
head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819
tW_per_site / (4 * 2.04e-9)
## 217,483.3
## mean: 212173

# par pre

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819

tW_per_site / (4 * 2.04e-9)
## 169,284.5

# par post

head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819
tW_per_site / (4 * 2.04e-9)
## 178,145.6
## mean: 173715
## (169284.5 + 178145.6)/2
