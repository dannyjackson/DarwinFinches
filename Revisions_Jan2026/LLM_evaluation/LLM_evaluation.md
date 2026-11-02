# LLM evaluation of GO terms (Figure 4)

Scripts and prompts to reproduce the LLM a priori relevance scores of GO terms and Figure 4.

## Folder layout
```
${PROJDIR}/LLM_evaluation/
  README.md
  prompts/
    all_pathways.txt          the list of GO terms to evaluate
    apriori_GO_Ranking_Prompt.txt   the prompt that ran 30 raters as sub-agents in one new conversation per source
  scripts/
    make_term_list.py         1. list the 1,008 GO terms tested in all three species (runs/all_pathways.txt)
    setup_runs.py             2. make the run folders: shuffled term list and filled rater prompt per run
    validate_runs.py          3. check every run's scores.tsv
    compile_matrix.py         4. combine a source's runs into one score matrix
    make_figure4.py           5. Figure 4, its source data and the top-terms table
  data/                       (inputs, not included)
    gsea/cra.fst.gsea.tsv  for.fst.gsea.tsv  par.fst.gsea.tsv
    score_matrices/score_matrix.<source>.tsv   (written by compile_matrix.py)
  runs/                       (written by the scripts and the raters)
  results/                    (written by make_figure4.py)
```
All paths default to `${PROJDIR}/LLM_evaluation/...` and can be changed with each script's options (`--help`).
If `PROJDIR` is not set, the scripts use the folder above `scripts/`.

## Requirements
Python 3 with numpy, pandas and matplotlib, plus Pillow for the TIFF. Font: Arial or Helvetica if installed,
otherwise Liberation Sans.

## Steps
```
export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/

# 1. term list given to the raters (GO ID and term name only)
python3 $SCRIPTDIR/make_term_list.py

# 2-4. for each source (one model on one account), e.g. Claude Opus 5.5 on the first author's account, batch a
SRC=informed.Opus5_5.run_a
python3 $SCRIPTDIR/setup_runs.py --source $SRC --n-runs 30
#    give each runs/$SRC/run_NNN/prompt.txt to a fresh, independent rater (see "Generating runs" below)
python3 $SCRIPTDIR/validate_runs.py --source $SRC        # re-run any failed run from scratch with a fresh rater
python3 $SCRIPTDIR/compile_matrix.py --source $SRC

# 5. Figure 4 from all score matrices
python3 $SCRIPTDIR/make_figure4.py
```

## Source names
Each source is named `<informed|naive>.[Provider.]<model>.<batch|account holder>`:
- **informed** sources ran on the first author's accounts. The last part is the batch (`run_a` or `run_b`), shown as Ia and Ib.
- **naive** sources ran on accounts of people unfamiliar with the study. The last part is the account holder's name, which the figure replaces with a code (`--codes Lexi=N1 Zach=N2 Rebecca=N3`; an unlisted holder gets the next free N number).
- Models: `Opus5_5` (Claude Opus 5.5), `Fable5_1` (Claude Fable 5.1), `gpt6_1sol` (GPT 6.1 Sol), `gpt6astra` (GPT 6 Astra). The optional provider part (e.g. `ChatGPT`) is ignored.

Examples: `informed.Opus5_5.run_a`, `naive.Fable5_1.Rebecca`, `naive.ChatGPT.gpt6astra.Rebecca`.

The figure uses 13 sources of 30 runs each (390 runs). Each model ran in two informed batches (Ia and Ib) and on naive accounts: Opus 5.5 on N1 and N2, and Fable 5.1, GPT 6.1 Sol and GPT 6 Astra on N3.

## Generating runs
Each run is one independent rating of all 1,008 terms by a rater that sees only `prompt.txt` and its own
`input_terms.txt`. Term order is shuffled per run with `random.Random(1000 + run)`, and the seeds are logged in
`runs/<source>/run_manifest.tsv`.

For Claude, each source was run in one new conversation with no project attached: it received
`prompts/coordinating_prompt.txt` and `all_pathways.txt`, and launched one sub-agent per run in batches of five.
STEPS 1–3 of that prompt do the same thing as `setup_runs.py`, `validate_runs.py` and `compile_matrix.py`. Rename
the matrix it returns to `score_matrix.<source>.tsv` and put it in `data/score_matrices/`.

[If the ChatGPT runs used a different procedure or coordinating prompt, describe it here.]

## Figure 4 (`make_figure4.py`)
- **A** Mean quadratic-weighted Cohen's κ between runs of each pair of sources. The diagonal is the mean over pairs of different runs within one source.
- **B** Distribution of scores for two focal terms (`--focal GO:0007596 GO:0002250`) in each source, with the mean score and the term's rank among all 1,008 terms by mean score.
- **C** Consensus terms: median score ≥4 across all runs pooled (`--consensus-score 4`). The panel shows those whose NES averaged across the three species is in the top 10% of all tested terms (`--nes-top 0.10`). For each term it shows the NES in each species, the mean NES, the rank by mean NES, the mean score and the gene-set size.

Species are labelled VF (vegetarian finch, `cra`), MGF (medium ground finch, `for`) and STF (small tree finch,
`par`). Change the names with `--names` and the fgsea tables with `--gsea VF=path MGF=path STF=path`.

Outputs in `results/`:
- `Fig4_LLM_consistency.pdf` (vector, fonts embedded), `.png` and `.tiff` (600 dpi). The figure is 18.3 cm wide (Science full width).
- `Fig4_LLM_consistency_caption.md` and `.txt`: caption with the current counts.
- `panelA_kappa_matrix.tsv`, `panelB_focal_terms.tsv` and `panelC_consensus_terms_NES.tsv`: the data behind each panel. The last lists all consensus terms with their NES, ranks, scores and fgsea padj; `high_NES` marks the terms shown.
- `table_top_terms.tsv`: the 10 terms with the highest mean NES (`--top-n 10`), with their rank in each species, median and mean score, percentage of runs scoring ≥4, and whether each is a consensus term. This is the table behind the results text.
- `sources.tsv` and `account_codes.tsv`: the sources and the account-code key.
