# LLM evaluation of GO terms (Figure 4)

Scripts and prompts to reproduce the LLM a priori relevance scores of GO terms and Figure 4.

## Folder layout
Inputs live in this repository (`${SCRIPTDIR}`). Results are written to `${PROJDIR}/LLM_evaluation/results/`,
outside the repository. Run everything from `${PROJDIR}`.
```
${SCRIPTDIR} = ${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
  README.md
  prompts/
    apriori_GO_Ranking_Prompt.txt   prompt given to each source's coordinating conversation
    all_pathways.txt                the 1,008 GO terms given to the raters (GO ID and term name)
  data/
    cra.fst.gsea.tsv  for.fst.gsea.tsv  par.fst.gsea.tsv   fgsea tables
    matrices/score_matrix.<source>.tsv                      one per source (30 runs each)
    accounts.tsv                                            which account each source ran on (for panel B)
  scripts/
    make_term_list.py   1. check that all_pathways.txt matches the fgsea tables (writes it only if missing)
    make_figure4.py     2. Figure 4, its source data and the top-terms table
    plot_top_NES_LLM.py 3. the top-NES plots from D1_run_fgsea.R, with the terms all models scored ≥4 in colour

${PROJDIR}/LLM_evaluation/results/   written by make_figure4.py
```
Every path can be changed with each script's options (`--help`).

## Steps
```
export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
export SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
cd $SCRIPTDIR
git pull
cd $PROJDIR

module load micromamba
# micromamba create -y -n llm_eval -c conda-forge python=3.12 numpy pandas scipy matplotlib pillow && micromamba activate llm_eval
micromamba activate llm_eval

# Make Figure 4 from all score matrices in data/matrices/

python3 $SCRIPTDIR/scripts/make_figure4.py

# 3. <species>.fst.top_NES.LLM.pdf: top 20 terms by NES, terms scored >=4 by every kept model coloured
python3 $SCRIPTDIR/scripts/plot_top_NES_LLM.py
```