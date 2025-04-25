θ = 2N e(f)μ, where μ is the substitution rate per generation with estimations of thet

# From Lamichhaney et al. 2015:
# We used the following previously reported estimated mutation rates for nuclear and mtDNA: nuclear DNA, 2.04 × 10−9 per site per year estimated from the synonymous mutation rate on the Darwin’s finches’ lineage since the split from zebra finch45; mtDNA, a fossil-calibrated divergence rate of 2.1% per million years for bird cytochrome b sequences46.

# Compute average Watterson's theta for each pop
# cra pre
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapre/crapre.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 

# cra post
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/crapost/crapost.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 

# for pre
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpre/forpre.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 

# for post
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/forpost/forpost.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 

# par pre
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpre/parpre.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 

# par post
grep 'NC_' /xdisk/mcnew/finches/dannyjackson/finches/analyses/thetas/parpost/parpost.thetas.idx.pestPG | grep -v 'NC_044601' | awk -v mu=2.04e-9 'NR > 1 { tw += $4; sites += $14 } END {
  tw_per_site = tw / sites;
  Ne = tw_per_site / (4 * mu);
  print "Genome-wide tW per site:", tw_per_site;
  print "Estimated Ne:", Ne
}' 