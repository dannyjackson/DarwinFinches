# Power of the empirical GSEA to detect parallel selection on blood coagulation

Question: if blood coagulation (GO:0007596) were under selection in all three species,
how often would our empirical GSEA pipeline rank it near the top?

## Design

For one simulated condition (Ne, generations, decline, s, f0), each iteration:

1. For each species (CRA, FOR, PAR), independently:
   - Start from that species' empirical gene list (`gene_scores/<sp>.fst.gene_scores.tsv`).
   - Give every gene an FST drawn from the neutral simulations with the same Ne, generations and decline.
     The neutral pool (~10,000 windows) is smaller than the gene list, so draws are with replacement.
   - Give every blood-coagulation gene in that species' ranking an FST drawn from the selected simulations
     (all genes in the term are under selection; drawn without replacement from the 1,000 selected windows).
   - Run fgsea exactly as `D1_run_fgsea.R`: `go_terms.tsv`, `scoreType = "pos"`, `minSize = 10`, `maxSize = 500`, `eps = 0`.
2. Combine species exactly as `D1_cross_species_analysis.R`: mean NES across the terms tested in all three species.
3. Record blood coagulation's rank and percentile (rank / number of terms × 100; smaller = stronger):
   per species (by NES and by nominal P, plus padj) and across species (by mean NES).

100 iterations per condition. A null control (s = 0: blood-coagulation genes also get neutral FST)
gives the false-positive baseline; its percentile should be roughly uniform.

Scripts (all in `Simulations_Oct2026/`):

- `gsea_power_empirical.R`: one condition, 100 iterations, resumable
- `run_gsea_power_empirical.sh`: SLURM wrapper for one condition
- `summarize_gsea_power_empirical.R`: summary table and power curves

## Set up

```
PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/Simulations_Oct2026
cd ${SCRIPTDIR} && git pull
cd ${PROJDIR}
mkdir -p ${PROJDIR}/output/logs
```

## Test and time one condition first (interactive)

```
module load R
cd ${PROJDIR}

Rscript ${SCRIPTDIR}/gsea_power_empirical.R \
  Ne=100k Gen=5 decline=1.00 s=0.1 f0=0.05 \
  n_iter=2 nproc=4 out_dir=${PROJDIR}/gsea_empirical_power_test
```

Check the log for:
- the number of blood-coagulation genes found in each species' ranking (expect ~25)
- the empirical `mean_score` range vs. the percentage of negative simulated FST
- seconds per iteration (3 fgsea runs); set `--time` in the wrapper from this

## Run

Main set: Ne = 100k, no decline, both generation intervals, all s × f0, plus a null control per interval.

```
##########################################
NE=100k
##########################################
D=1.00
GEN_TIMES=(5 10)
SELS=(0.01 0.02 0.05 0.1 0.2 0.5)
F0S=(0.05 0.1 0.2 0.4 0.6)

for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done

D=0.80
for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done

##########################################
NE=50k
##########################################
D=1.00
GEN_TIMES=(5 10)
SELS=(0.01 0.02 0.05 0.1 0.2 0.5)
F0S=(0.05 0.1 0.2 0.4 0.6)

for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done


D=0.80

for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done

##########################################
NE=200k
##########################################
D=1.00
GEN_TIMES=(5 10)
SELS=(0.01 0.02 0.05 0.1 0.2 0.5)
F0S=(0.05 0.1 0.2 0.4 0.6)

for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done


D=0.80

for G in "${GEN_TIMES[@]}"; do
  # null control
  sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s 0.0

  for S in "${SELS[@]}"; do
    for F0 in "${F0S[@]}"; do
      sbatch ${SCRIPTDIR}/run_gsea_power_empirical.sh -N ${NE} -o ${G} -r ${D} -s ${S} -f ${F0}
      sleep 0.1
    done
  done
done
```

To extend to the other demographies, loop over `NE` in `50k 100k 200k` and `D` in `1.00 0.8`
(each Ne × decline × Gen needs its own null control).

Jobs are resumable: if one times out, resubmit it with the same flags.

## Summarise and plot

```
cd ${PROJDIR}
module load R
Rscript ${SCRIPTDIR}/summarize_gsea_power_empirical.R
```

Outputs in `gsea_empirical_power/`:

- `summary.tsv`, one row per condition:
  - `median_pct_cross`: median percentile of blood coagulation by mean NES across species
  - `power_cross_top1` / `_top5` / `_top10`: proportion of iterations in the top 1/5/10% by mean NES
  - `power_sig_all3` / `power_sig_any`: proportion with padj < 0.05 in all three / at least one species
  - `prop_positive_all3`: proportion with NES > 0 in all three species
  - `power_sig_<SP>`, `median_pct_<SP>`: per-species results
- `plots/<metric>.Gen_<g>.png`: power curves (x = s, facets = f0); dashed line = the s = 0 control
