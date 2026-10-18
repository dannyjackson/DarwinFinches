# B2. Effective population size from Watterson's θ

Ne is estimated for each population (3 species × pre/post invasion) from genome-wide autosomal Watterson's θ:

```
θ_W (per site) = Σ tW / Σ nSites                (autosomes only)
μ_gen          = μ_year × generation time
Ne             = θ_W / (4 μ_gen)
```

**Mutation rate.** From Lamichhaney et al. 2015: *"nuclear DNA, 2.04 × 10⁻⁹ per site per year estimated from the synonymous mutation rate on the Darwin's finches' lineage since the split from zebra finch."* This is a rate **per year**. θ = 4Neμ needs μ **per generation**, so the yearly rate is multiplied by the generation time: 5 years gives μ_gen = 1.02 × 10⁻⁸. The MSMC plots are scaled with μ = 1 × 10⁻⁸ and g = 4.5, which is consistent.

**What changed from the previous version of this file**

| Previous | Now |
|---|---|
| `Ne = θ / (4 × 2.04e-9)`: the yearly rate used as per-generation, so the result is Ne × generation time (≈5× too high) | `Ne = θ / (4 × 2.04e-9 × g)`, with g = 5 as the primary value plus a sensitivity table (3, 4.5, 7 yr) |
| crapre divided by `sum(nSites)`; the other five by a hard-coded `946505819` | every population divides by its own `sum(nSites)` |
| crapre θ came from `datafiles/safs` (the SNP-only SAFs used for FST) | all six populations come from `analyses/effective_pop/safs` (all sites, autosomes; B2_1) |
| crapost dropped from the R loop | all six populations processed |
| "Results" block listed the pre/post means as the FOR and PAR "post" values (the source of the wrong 171,500 and 209,518 in the manuscript) | species means are computed in the script from the per-population values |
| no uncertainty | weighted delete-one-chromosome jackknife SE (Busing et al. 1999) |

No new environment is needed: step 3 uses only the Python 3 standard library.

---

## 1. Site allele frequency likelihoods (all sites, autosomes)

`B2_1_run_thetas_array.sh` (unchanged): `angsd -dosaf 1 -GL 1 -rf autosomes.txt -anc ${REF}` with `-setMinDepthInd 4 -minInd 4 -minQ 30 -minMapQ 30` and **no `-sites` or MAF filter**. Invariant sites must be included for θ per site to be meaningful.

```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop
mkdir -p slurm_output
sbatch --array=0-5 ~/programs/DarwinFinches/Scripts/B_PopulationStructureAnalyses/B2_1_run_thetas_array.sh
```

Array index → population: 0 crapre, 1 crapost, 2 forpre, 3 forpost, 4 parpre, 5 parpost.

## 2. SFS, per-site θ, per-chromosome summaries

`B2_2_run_thetas_array.sh` (corrected: all six populations, `-P ${THREADS}`): `realSFS` → `realSFS saf2theta` → `thetaStat do_stat`. It writes `thetas/${POP}.thetas.idx.pestPG` with one row per chromosome.

```
sbatch --array=0-5 ~/programs/DarwinFinches/Scripts/B_PopulationStructureAnalyses/B2_2_run_thetas_array.sh
```

Notes:
- The SFS is unfolded with the reference as ancestral (`-anc ${REF}`). Watterson's θ depends only on whether a site is segregating, not on which allele is ancestral, so tW is unaffected. Folding would matter for Fay & Wu's H, not for this.
- Do **not** add `-win`/`-step` to `do_stat` here. Step 3 expects one row per chromosome and stops if it finds duplicate chromosome names.

Before step 3, check that every population finished and has a row per autosome:

```
for POP in crapre crapost forpre forpost parpre parpost; do
  printf "%s\t" $POP; grep -c '^(' thetas/${POP}.thetas.idx.pestPG
done
```

## 3. Ne with a per-generation mutation rate

```
cd /xdisk/mcnew/finches/dannyjackson/finches/analyses/effective_pop

python3 ~/programs/DarwinFinches/Scripts/B_PopulationStructureAnalyses/B2_3_compute_Ne.py \
  --theta_dir thetas \
  --mu_year 2.04e-9 \
  --gen_time 5 \
  --gen_sens 3 4.5 7 \
  --exclude NC_044601 \
  --prefix NC_ \
  --out Ne_estimates.tsv
```

(If `python3` is not on the path on the login/compute node: `module load python`.)

What it does:
- Keeps chromosomes whose names start with `NC_` (drops unplaced scaffolds) and drops the Z (`NC_044601`, with or without a version suffix). The SAFs were already restricted with `-rf autosomes.txt`, so this is a safeguard.
- θ_W per site = Σ tW / Σ nSites over the kept chromosomes, separately for each population.
- Ne = θ_W / (4 × μ_year × g) for g = 5 (primary) and 3, 4.5, 7 (sensitivity).
- SE: weighted delete-one-chromosome jackknife (Busing et al. 1999), which accounts for unequal chromosome lengths.
- Species values: mean of pre and post.
- Also writes `Ne_legacy_mu_per_year` (the old calculation) so the change can be checked: legacy / Ne_g5 = 5 exactly.

Output columns in `Ne_estimates.tsv`: `pop species timepoint n_chrom n_sites theta_W_per_site theta_W_SE Ne_g5 Ne_g5_SE Ne_g3 … Ne_g4.5 … Ne_g7 … Ne_legacy_mu_per_year`.

## 4. Results

*(Fill in from `Ne_estimates.tsv` after running steps 1–3.)*

| Species | Pre Ne (± SE) | Post Ne (± SE) | Mean |
|---|---|---|---|
| Vegetarian finch (CRA) | | | |
| Medium ground finch (FOR) | | | |
| Small tree finch (PAR) | | | |

For orientation only, the previous values divided by 5 are:
- CRA ≈ 31,200 (pre) / 18,100 (post)
- FOR ≈ 41,400 / 43,500
- PAR ≈ 33,900 / 35,600

The rerun will differ for CRA in particular, because crapre previously came from a different SAF set and denominator.

## 5. Text for the manuscript (update numbers after rerunning)

Supplement:

> Effective population size was estimated from Watterson's θ as Ne = θ_W/(4μ), where θ_W is the autosomal mean per-site Watterson's θ estimated in ANGSD (realSFS saf2theta and thetaStat, using all sites passing filters) and μ is the per-generation mutation rate. We converted the substitution rate of 2.04 × 10⁻⁹ per site per year (Lamichhaney et al. 2015) to a per-generation rate assuming a generation time of 5 years (μ = 1.02 × 10⁻⁸ per site per generation). Standard errors were estimated with a weighted delete-one-chromosome jackknife. Estimates under alternative generation times (3–7 years) are given in Table S[x].

Main text: replace "123,340, 171,500, and 209,518" with the species means from `Ne_estimates.tsv` (Ne_g5).

The simulation grid (Ne = 50,000, 100,000, 200,000) and the supplement sentence "These estimates produced Ne values of approximately 90,000–212,000" both need revisiting once these numbers are in. The corrected Ne is expected to fall at about 18,000–44,000.
