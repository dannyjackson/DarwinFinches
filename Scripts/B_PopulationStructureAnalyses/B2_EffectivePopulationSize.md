θ = 2N e(f)μ, where μ is the substitution rate per generation with estimations of thet

# From Lamichhaney et al. 2015:
# We used the following previously reported estimated mutation rates for nuclear and mtDNA: nuclear DNA, 2.04 × 10−9 per site per year estimated from the synonymous mutation rate on the Darwin’s finches’ lineage since the split from zebra finch45; mtDNA, a fossil-calibrated divergence rate of 2.1% per million years for bird cytochrome b sequences46.

# Compute average Watterson's theta for each pop
# cra pre
head -n 1 /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG > /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG.autosomes
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG | grep -v 'NC_044601' >> /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG.autosomes

R

df <- read.csv("/xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG.autosomes", sep='\t')
tW_per_site <- sum(df$tW) / 946505819

tW_per_site / (4 * 2.04e-9)

## 99,978.85

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
