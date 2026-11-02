# LLM evaluation of GO terms (Figure 4)

Scripts and prompts to reproduce the LLM a priori relevance scores of GO terms and Figure 4.

## Folder layout
Inputs live in this repository (`${SCRIPTDIR}`). Runs and results are written to `${PROJDIR}/LLM_evaluation/`,
outside the repository. Run everything from `${PROJDIR}`.
```
${SCRIPTDIR} = ${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
  README.md
  prompts/
    apriori_GO_Ranking_Prompt.txt   coordinating prompt; the rater prompt is the text between its two '---' lines
    all_pathways.txt                the 1,008 GO terms given to the raters (GO ID and term name)
  data/
    cra.fst.gsea.tsv  for.fst.gsea.tsv  par.fst.gsea.tsv   fgsea tables
    matrices/score_matrix.<source>.tsv                      one per source (30 runs each)
  scripts/
    make_term_list.py   1. build all_pathways.txt from the fgsea tables, or check the existing one matches
    setup_runs.py       2. make the run folders: shuffled term list and filled rater prompt per run
    validate_runs.py    3. check every run's scores.tsv
    compile_matrix.py   4. combine a source's runs into data/matrices/score_matrix.<source>.tsv
    make_figure4.py     5. Figure 4, its source data and the top-terms table

${PROJDIR}/LLM_evaluation/
  runs/<source>/run_NNN/   written by setup_runs.py and the raters
  results/                 written by make_figure4.py
```
The scripts only write inside the repository in two cases. `make_term_list.py` writes `prompts/all_pathways.txt`
only if that file is missing (otherwise it checks that the file matches the fgsea tables). `compile_matrix.py`
writes new score matrices to `data/matrices/` and will not overwrite one without `--force`. Every path can be changed
with each script's options (`--help`).

## Requirements
Python 3 with numpy, pandas and matplotlib, plus Pillow for the TIFF. Font: Arial or Helvetica if installed,
otherwise Liberation Sans.

## Steps
```
export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
export SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
cd $SCRIPTDIR
git pull
cd $PROJDIR

# 1. term list given to the raters (GO ID and term name only)
python3 $SCRIPTDIR/scripts/make_term_list.py

# 2. Figure 4 from all score matrices in data/matrices/
python3 $SCRIPTDIR/scripts/make_figure4.py
```

## Source names
Each source is named `<informed|naive>.[Provider.]<model>.<batch|account holder>`:
- **informed** sources ran on the first author's accounts. The last part is the batch (`run_a` or `run_b`), shown as Ia and Ib.
- **naive** sources ran on accounts of people unfamiliar with the study. The last part is the account holder's name, which the figure replaces with a code (`--codes Lexi=N1 Zach=N2 Rebecca=N3 Nora=N4`, the default; an unlisted holder gets the next free N number).
- Models: `Opus5_5` (Claude Opus 5.5), `Fable5_1` (Claude Fable 5.1), `gpt6_1sol` (GPT 6.1 Sol), `gpt6astra` (GPT 6 Astra). The optional provider part (e.g. `ChatGPT`) is ignored.

Examples: `informed.Opus5_5.run_a`, `naive.Fable5_1.Rebecca`, `naive.ChatGPT.gpt6astra.Rebecca`.

Each model ran in two informed batches (Ia and Ib) and on naive accounts: Opus 5.5 on N1 and N2; Fable 5.1 and GPT 6.1 Sol on N3 and N4; GPT 6 Astra on N3. `make_figure4.py` prints the number of sources and runs it read, and skips (with a warning) any matrix that is byte-identical to another.

## Generating runs
Each run is one independent rating of all 1,008 terms by a rater that sees only `prompt.txt` and its own
`input_terms.txt`. Term order is shuffled per run with `random.Random(1000 + run)`, and the seeds are logged in
`runs/<source>/run_manifest.tsv`.

For Claude, each source was run in one new conversation with no project attached: it received
`prompts/apriori_GO_Ranking_Prompt.txt` and `prompts/all_pathways.txt`, and launched one sub-agent per run in batches
of five. STEPS 1–3 of that prompt do the same thing as `setup_runs.py`, `validate_runs.py` and `compile_matrix.py`.
Rename the matrix it returns to `score_matrix.<source>.tsv` and put it in `data/matrices/`.

[If the ChatGPT runs used a different procedure or coordinating prompt, describe it here.]

## Figure 4 (`make_figure4.py`)
- **A** Mean quadratic-weighted Cohen's κ between runs of each pair of sources. The diagonal is the mean over pairs of different runs within one source.
- **B** Distribution of scores for two focal terms (`--focal GO:0007596 GO:0002250`) in each source, with the mean score and the term's rank among all 1,008 terms by mean score.
- **C** Consensus terms: median score ≥4 across all runs pooled (`--consensus-score 4`). The panel shows those whose NES averaged across the three species is in the top 10% of all tested terms (`--nes-top 0.10`). For each term it shows the NES in each species, the mean NES, the rank by mean NES, the mean score and the gene-set size.

Species are labelled VF (vegetarian finch, `cra`), MGF (medium ground finch, `for`) and STF (small tree finch,
`par`). Change the names with `--names` and the fgsea tables with `--gsea VF=path MGF=path STF=path`.

Outputs in `${PROJDIR}/LLM_evaluation/results/`:
- `Fig4_LLM_consistency.pdf` (vector, fonts embedded), `.png` and `.tiff` (600 dpi). The figure is 18.3 cm wide (Science full width).
- `Fig4_LLM_consistency_caption.md` and `.txt`: caption with the current counts.
- `panelA_kappa_matrix.tsv`, `panelB_focal_terms.tsv` and `panelC_consensus_terms_NES.tsv`: the data behind each panel. The last lists all consensus terms with their NES, ranks, scores and fgsea padj; `high_NES` marks the terms shown.
- `table_top_terms.tsv`: the 10 terms with the highest mean NES (`--top-n 10`), with their rank in each species, median and mean score, percentage of runs scoring ≥4, and whether each is a consensus term. This is the table behind the results text.
- `sources.tsv` and `account_codes.tsv`: the sources and the account-code key.
