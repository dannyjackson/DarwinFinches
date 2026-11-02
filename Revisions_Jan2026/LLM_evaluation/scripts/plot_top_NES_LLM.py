#!/usr/bin/env python3
"""
Remake the <species>.fst.top_NES.pdf plots from D1_run_fgsea.R (the 20 GO terms with the highest NES; point size =
gene-set size, triangle = nominal P < 0.05), colouring the terms that every kept LLM model scored >= 4 (median score
across that model's runs; the centre of Fig. 4E).

Run make_figure4.py first: this script reads its table_venn_terms.tsv.

Inputs (defaults):
  ${SCRIPTDIR}/data/{cra,for,par}.fst.gsea.tsv            fgsea tables (pathway, term_name, NES, pval, size)
  ${PROJDIR}/LLM_evaluation/results/table_venn_terms.tsv   from make_figure4.py (column all_models_ge4)

Usage:
  export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
  export SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
  cd $PROJDIR
  python3 $SCRIPTDIR/scripts/make_figure4.py
  python3 $SCRIPTDIR/scripts/plot_top_NES_LLM.py

Outputs (in --out, default ${PROJDIR}/LLM_evaluation/results):
  <species>.fst.top_NES.LLM.pdf and .png for each fgsea table
"""
import argparse, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CODE = os.environ.get("SCRIPTDIR") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(CODE, "data")
WORK = os.path.join(os.environ.get("PROJDIR", os.getcwd()), "LLM_evaluation")
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--gsea", nargs="+", default=[os.path.join(DATA, f + ".fst.gsea.tsv") for f in ("cra", "for", "par")])
ap.add_argument("--venn", default=os.path.join(WORK, "results", "table_venn_terms.tsv"),
                help="table_venn_terms.tsv from make_figure4.py")
ap.add_argument("--out", default=os.path.join(WORK, "results"))
ap.add_argument("--top", type=int, default=20, help="number of terms with the highest NES to plot")
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)

if not os.path.exists(A.venn): sys.exit(f"{A.venn} not found; run make_figure4.py first")
V = pd.read_csv(A.venn, sep="\t")
models = sorted({m.strip() for s in V.models_ge4 for m in s.split(",") if m.strip() != "none"})
SHARED = set(V.loc[V.all_models_ge4.astype(str).str.lower() == "true", "pathway"])
print(f"{len(SHARED)} terms scored >= 4 by every model ({', '.join(models)})")

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
                     "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
HI, LO = "#d0743c", "#000000"          # coloured = scored >= 4 by every model
for f in A.gsea:
    g = pd.read_csv(f, sep="\t").dropna(subset=["NES", "pval"])
    g["term_name"] = g.term_name.fillna(g.pathway)
    d = g.sort_values("NES", ascending=False).head(A.top).reset_index(drop=True)
    d["shared"] = d.pathway.isin(SHARED)
    y = np.arange(len(d))[::-1]                      # highest NES at the top, as in D1_run_fgsea.R
    smin, smax = g["size"].min(), g["size"].max()
    area = lambda s: 20 + 180 * (s - smin) / max(smax - smin, 1)

    fig, ax = plt.subplots(figsize=(8, 7))
    for sig, mk in ((False, "o"), (True, "^")):
        for sh, col in ((False, LO), (True, HI)):
            k = (d.pval < 0.05).eq(sig) & d.shared.eq(sh)
            ax.scatter(d.NES[k], y[k], s=area(d["size"][k]), marker=mk, color=col, zorder=3)
    ax.set_yticks(y, d.term_name)
    for lab, sh in zip(ax.get_yticklabels(), d.shared):
        if sh: lab.set_color(HI); lab.set_fontweight("bold")
    ax.set_xlabel("Normalized enrichment score (NES)")
    ax.grid(color="#e5e5e5", lw=.6); ax.set_axisbelow(True)
    ax.set_ylim(-0.7, len(d) - 0.3)

    sizes = [s for s in (10, 50, 100, 200, 400) if smin <= s <= smax] or [int(smin), int(smax)]
    h_col = [Line2D([], [], marker="o", ls="", color=HI, ms=7, label=f"Scored ≥4 by all\n{len(models)} models"),
             Line2D([], [], marker="o", ls="", color=LO, ms=7, label="Other terms")]
    h_shp = [Line2D([], [], marker="o", ls="", color="#555555", ms=7, label="P ≥ 0.05"),
             Line2D([], [], marker="^", ls="", color="#555555", ms=7, label="P < 0.05")]
    h_size = [plt.scatter([], [], s=area(s), color="#555555", label=str(s)) for s in sizes]
    l1 = ax.legend(handles=h_col, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, title="A priori relevance")
    ax.add_artist(l1)
    l2 = ax.legend(handles=h_shp, loc="upper left", bbox_to_anchor=(1.02, 0.78), frameon=False, title="Nominal P")
    ax.add_artist(l2)
    l3 = ax.legend(handles=h_size, loc="upper left", bbox_to_anchor=(1.02, 0.58), frameon=False, title="Gene-set size",
              labelspacing=1.1)
    for l in (l1, l2, l3): l._legend_box.align = "left"
    fig.tight_layout()
    base = os.path.join(A.out, os.path.basename(f).replace(".gsea.tsv", "") + ".top_NES.LLM")
    kw = dict(bbox_inches="tight", bbox_extra_artists=[l1, l2, l3])   # include all three legends in the saved area
    fig.savefig(base + ".pdf", **kw); fig.savefig(base + ".png", dpi=300, **kw); plt.close(fig)
    print(f"{os.path.basename(f)}: {int(d.shared.sum())} of the top {len(d)} terms by NES are shared -> {base}.pdf")
