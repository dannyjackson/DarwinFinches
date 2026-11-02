#!/usr/bin/env python3
"""
Figure 4: consistency of LLM a priori GO-term relevance scores, and the consensus terms with the highest
enrichment.

  A  Mean run-vs-run quadratic-weighted kappa between every pair of sources (diagonal = within source)
  B  Score distributions for two focal terms (default: blood coagulation, adaptive immune response)
  C  Consensus terms (median score >= 4 across all runs pooled) whose NES averaged across species is in the
     top 10% of all tested terms
Also writes table_top_terms.tsv: the --top-n terms with the highest mean NES and their LLM scores (results text).

Inputs (defaults are relative to ${PROJDIR}/LLM_evaluation; every path can be overridden):
  data/score_matrices/score_matrix.<informed|naive>.[Provider.]<model>.<batch|holder>.tsv
        e.g. score_matrix.informed.Opus5_5.run_a.tsv, score_matrix.naive.ChatGPT.gpt6_1sol.Rebecca.tsv
        columns: pathway, term_name, run_001 ... run_N  (scores 1-5), as written by compile_matrix.py
        informed batches run_a / run_b (or run2) are coded Ia / Ib; naive holders are coded via --codes
  data/gsea/{cra,for,par}.fst.gsea.tsv  fgsea tables (columns pathway, NES, padj, size, ...),
        shown as VF (vegetarian finch), MGF (medium ground finch), STF (small tree finch)

Usage:
  export PROJDIR=/path/to/project
  python3 ${PROJDIR}/LLM_evaluation/scripts/make_figure4.py
  python3 ${PROJDIR}/LLM_evaluation/scripts/make_figure4.py --codes Lexi=N1 Zach=N2 Rebecca=N3 NewPerson=N4

Outputs (in --out, default ${PROJDIR}/LLM_evaluation/results):
  Fig4_LLM_consistency.pdf (vector, fonts embedded), .png and .tiff (600 dpi), caption (.md, .txt)
  panelA_kappa_matrix.tsv, panelB_focal_terms.tsv, panelC_consensus_terms_NES.tsv, table_top_terms.tsv,
  sources.tsv, account_codes.tsv

Requires: Python 3, numpy, pandas, matplotlib (Pillow for TIFF).
"""
import argparse, glob, hashlib, os, re, sys, warnings
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
warnings.filterwarnings("ignore")

ROOT = os.path.join(os.environ["PROJDIR"], "LLM_evaluation") if "PROJDIR" in os.environ else \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--matrices", default=os.path.join(ROOT, "data", "score_matrices"))
ap.add_argument("--gsea", nargs="+", default=[f"{lab}={os.path.join(ROOT, 'data', 'gsea', f + '.fst.gsea.tsv')}"
                                               for lab, f in (("VF", "cra"), ("MGF", "for"), ("STF", "par"))],
                help="LABEL=path for each species; the label is what the figure shows")
ap.add_argument("--out", default=os.path.join(ROOT, "results"))
ap.add_argument("--codes", nargs="*", default=["Lexi=N1", "Zach=N2", "Rebecca=N3"],
                help="codes for naive account holders (named in the file name); unlisted ones get the next free N number")
ap.add_argument("--focal", nargs="+", default=["GO:0007596", "GO:0002250"])
ap.add_argument("--names", nargs="*", default=[], help="LABEL=full name overrides for the species legend, e.g. VF='Vegetarian finch'")
ap.add_argument("--consensus-score", type=int, default=4, help="main panel C: score threshold for consensus terms (all runs pooled)")
ap.add_argument("--nes-top", type=float, default=0.10, help="main panel C: show consensus terms whose mean NES is in this top fraction of all terms")
ap.add_argument("--top-n", type=int, default=10, help="number of top terms by mean NES listed in table_top_terms.tsv")
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)

# ------------------------------------------------------------------ load matrices
files = sorted(glob.glob(os.path.join(A.matrices, "score_matrix.*.tsv")))
if not files: sys.exit(f"No score_matrix.*.tsv files in {A.matrices}")
# File names: score_matrix.<informed|naive>.[Provider.]<model>.<batch or account holder>.tsv
#   informed: last token = run_a / run_b / run2 (batch a or b on the first author's accounts)
#   naive:    last token = the account holder's name (coded N1, N2 ... in the figure)
MODEL_NAMES = {"opus5_5": "Opus 5.5", "fable5_1": "Fable 5.1", "gpt6_1sol": "GPT 6.1 Sol", "gpt6astra": "GPT 6 Astra"}
MODEL_SHORT = {"opus5_5": "Opus", "fable5_1": "Fable", "gpt6_1sol": "Sol", "gpt6astra": "Astra"}
MODEL_ORDER = ["opus5_5", "fable5_1", "gpt6_1sol", "gpt6astra"]
def pretty_model(m): return MODEL_NAMES.get(m.lower(), m)
def short_model(m): return MODEL_SHORT.get(m.lower(), m)
def model_rank(m): return MODEL_ORDER.index(m.lower()) if m.lower() in MODEL_ORDER else 99
def batch_letter(tok):
    t = tok.lower()
    if t in ("run_a", "runa", "run1", "run_1", "a"): return "a"
    if t in ("run_b", "runb", "run2", "run_2", "b"): return "b"
    m = re.match(r"run_?([a-z])$", t)
    return m.group(1) if m else t

meta, mats, seen, ids, names = [], {}, {}, None, None
for f in files:
    h = hashlib.md5(open(f, "rb").read()).hexdigest()
    b = os.path.basename(f)[len("score_matrix."):-4]
    tok = b.split(".")
    if tok[0].lower() not in ("informed", "naive") or len(tok) < 3:
        sys.exit(f"{b}: expected score_matrix.<informed|naive>.[Provider.]<model>.<batch|name>.tsv")
    kind, model, last = tok[0].lower(), tok[-2], tok[-1]
    if h in seen:
        print(f"skip {b}: identical to {seen[h]}"); continue
    seen[h] = b
    d = pd.read_csv(f, sep="\t", dtype={"pathway": str}).set_index("pathway")
    if ids is None:
        ids = list(d.index); names = d["term_name"].to_dict() if "term_name" in d else {i: i for i in ids}
    runs = [c for c in d.columns if c.startswith("run_")]
    X = d.reindex(ids)[runs]
    if X.isna().any().any() or set(d.index) != set(ids):
        sys.exit(f"{b}: GO IDs differ from the first matrix or values are missing")
    if not X.stack().between(1, 5).all(): sys.exit(f"{b}: scores outside 1-5")
    naive = kind == "naive"
    meta.append(dict(source=b, model=model.lower(), naive=naive, holder=last if naive else "first author",
                     batch="" if naive else batch_letter(last), n_runs=len(runs)))
    mats[b] = X.to_numpy(float)
meta = pd.DataFrame(meta)
CODES = dict(c.split("=", 1) for c in A.codes)
for hld in meta[meta.naive].holder.unique():        # auto-code naive holders not given in --codes
    if hld in CODES: continue
    used = [int(v[1:]) for v in CODES.values() if v[1:].isdigit()]
    CODES[hld] = f"N{max(used, default=0) + 1}"
meta["code"] = [CODES[h] if n else f"I{bt}" for h, n, bt in zip(meta.holder, meta.naive, meta.batch)]
meta["label"] = [f"{pretty_model(m)} · {c}" for m, c in zip(meta.model, meta.code)]
meta["order"] = [(model_rank(m), n, c) for m, n, c in zip(meta.model, meta.naive, meta.code)]
meta = meta.sort_values("order").drop(columns="order").reset_index(drop=True)
pd.DataFrame({"account_holder": ["first author"] + list(CODES), "code": ["Ia / Ib (two batches)"] + list(CODES.values())}
             ).to_csv(os.path.join(A.out, "account_codes.tsv"), sep="\t", index=False)
SRC = meta.source.tolist(); LAB = dict(zip(meta.source, meta.label))
meta.to_csv(os.path.join(A.out, "sources.tsv"), sep="\t", index=False)
print(f"{len(SRC)} sources, {meta.n_runs.sum()} runs")

# ------------------------------------------------------------------ panel A: kappa matrix
def qwk(X, Y):
    n = X.shape[0]
    obs = (X ** 2).mean(0)[:, None] + (Y ** 2).mean(0)[None] - 2 * (X.T @ Y) / n
    exp = X.var(0)[:, None] + Y.var(0)[None] + (X.mean(0)[:, None] - Y.mean(0)[None]) ** 2
    return 1 - obs / exp
K = pd.DataFrame(index=SRC, columns=SRC, dtype=float)
for a in SRC:
    for b in SRC:
        k = qwk(mats[a], mats[b])
        K.loc[a, b] = k[np.triu_indices(k.shape[0], 1)].mean() if a == b else k.mean()
K.round(4).to_csv(os.path.join(A.out, "panelA_kappa_matrix.tsv"), sep="\t")

# ------------------------------------------------------------------ NES from the fgsea tables
gs = dict(g.split("=", 1) for g in A.gsea)
for sp, f in gs.items():
    if not os.path.exists(f): sys.exit(f"GSEA table for {sp} not found: {f} (set --gsea SP=path ...)")
G = {sp: pd.read_csv(f, sep="\t").set_index("pathway").reindex(ids) for sp, f in gs.items()}
NES = pd.DataFrame({sp: g.NES for sp, g in G.items()})
if NES.isna().any().any(): sys.exit("Some GO IDs are missing from a GSEA table")
NES["Mean"] = NES.mean(axis=1)
SPECIES = list(NES.columns)

ALLX = np.hstack([mats[s] for s in SRC])   # all runs pooled: terms x runs

# ------------------------------------------------------------------ panel B: focal terms
F = []
for s in SRC:
    X = mats[s]; mean = X.mean(1); rank = pd.Series(mean).rank(ascending=False, method="min").to_numpy()
    for go in A.focal:
        i = ids.index(go); x = X[i]
        F.append(dict(source=s, label=LAB[s], go=go, term=names.get(go, go), n_runs=len(x),
                      **{f"n_{k}": int((x == k).sum()) for k in range(1, 6)},
                      mean=x.mean(), rank_in_source=int(rank[i])))
F = pd.DataFrame(F); F.round(3).to_csv(os.path.join(A.out, "panelB_focal_terms.tsv"), sep="\t", index=False)


# ------------------------------------------------------------------ main panel C: consensus terms with high NES
# Consensus = median score of all runs pooled >= --consensus-score; high NES = mean NES across species in the
# top --nes-top fraction of all tested terms.
cons = np.median(ALLX, 1) >= A.consensus_score
mnes = NES["Mean"]
rank_nes = mnes.rank(ascending=False, method="min").astype(int)
cut = mnes.quantile(1 - A.nes_top)
size = G[SPECIES[0]]["size"] if "size" in G[SPECIES[0]] else pd.Series(np.nan, index=ids)
Cdf = pd.DataFrame({"pathway": ids, "term_name": [names[i] for i in ids], "size": size.to_numpy(),
                    "mean_score_all_runs": ALLX.mean(1), "median_score_all_runs": np.median(ALLX, 1),
                    **{f"NES_{sp}": NES[sp].to_numpy() for sp in SPECIES},
                    "rank_mean_NES": rank_nes.to_numpy(), "consensus": cons, "high_NES": (mnes >= cut).to_numpy()})
for sp in SPECIES[:-1]:
    if "padj" in G[sp]: Cdf[f"padj_{sp}"] = G[sp]["padj"].to_numpy()
Cdf = Cdf[Cdf.consensus].sort_values("rank_mean_NES")
Cdf.round(4).to_csv(os.path.join(A.out, "panelC_consensus_terms_NES.tsv"), sep="\t", index=False)
TOP = Cdf[Cdf.high_NES].reset_index(drop=True)
# table for the results text: top terms by mean NES and how the raters scored them
Tt = pd.DataFrame({"pathway": ids, "term_name": [names[i] for i in ids],
                   **{f"NES_{sp}": NES[sp].to_numpy() for sp in SPECIES},
                   "rank_mean_NES": rank_nes.to_numpy(),
                   **{f"rank_{sp}": NES[sp].rank(ascending=False, method="min").astype(int).to_numpy() for sp in SPECIES[:-1]},
                   "median_score": np.median(ALLX, 1), "mean_score": ALLX.mean(1),
                   "rank_mean_score": pd.Series(ALLX.mean(1)).rank(ascending=False, method="min").astype(int).to_numpy(),
                   "pct_runs_ge4": 100 * (ALLX >= 4).mean(1), "consensus": cons})
Tt = Tt.sort_values("rank_mean_NES").head(A.top_n)
Tt.round(3).to_csv(os.path.join(A.out, "table_top_terms.tsv"), sep="\t", index=False)
print(f"top {A.top_n} terms by mean NES: {int(Tt.consensus.sum())} consensus terms")
print(f"panel C: {len(Cdf)} consensus terms, {len(TOP)} with mean NES in the top {A.nes_top:.0%} (>= {cut:.2f})")

# ------------------------------------------------------------------ figure style (Science: 18.3 cm wide)
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 2.5, "pdf.fonttype": 42, "ps.fonttype": 42, "legend.fontsize": 6})
INK, MUTED, GRID = "#000000", "#4d4d4d", "#e0e0e0"
# species colours match the study-design figure (VF teal, MGF salmon, STF green); the cross-species mean is purple
SPC = dict(zip(SPECIES, ["#44afad", "#f37f7c", "#a5cb66", "#6250d6", "#eda100", "#e87ba4"]))
FULL = {"VF": "Vegetarian finch", "MGF": "Medium ground finch", "STF": "Small tree finch", "Mean": "Mean of species"}
FULL.update(dict(n.split("=", 1) for n in A.names))
def abbrev(t):
    t = t[0].upper() + t[1:]
    for a, b in (("Positive regulation", "Pos. reg."), ("Negative regulation", "Neg. reg."), ("Regulation", "Reg."),
                 ("positive regulation", "pos. reg."), ("negative regulation", "neg. reg.")):
        t = t.replace(a, b)
    return t
SPM = dict(zip(SPECIES, ["o", "s", "D", "^", "v", "P"]))
SEQ = ["#ededed", "#c6dbef", "#85b4de", "#3787c8", "#0b4a8c"]  # scores 1..5, one hue light->dark
W = 18.3 / 2.54

def letter_title(fig, x, y, letter, title):
    """Bold panel letter followed by a bold title, e.g. 'A. Run-vs-run agreement'."""
    rend = fig.canvas.get_renderer()
    t1 = fig.text(x, y, f"{letter}.", fontsize=8, fontweight="bold", va="bottom")
    w = t1.get_window_extent(rend).transformed(fig.transFigure.inverted()).width
    fig.text(x + w + 0.006, y, title, fontsize=8, fontweight="bold", va="bottom")

def save(fig, base):
    fig.savefig(base + ".pdf"); fig.savefig(base + ".png", dpi=600)
    try: fig.savefig(base + ".tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    except Exception as e: print("TIFF not written:", e)

# ================================================================== MAIN FIGURE: A, B, C
fig = plt.figure(figsize=(W, W * 0.86))
outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.36 + 0.045 * max(len(TOP), 1)], hspace=0.62,
                         left=0.105, right=0.935, top=0.95, bottom=0.085)
gsp = outer[0].subgridspec(1, 2, width_ratios=[1.25, 0.86], wspace=0.52)
model_of = meta.set_index("source").model.to_dict()
breaks = [i for i in range(1, len(SRC)) if model_of[SRC[i]] != model_of[SRC[i - 1]]]

# A heatmap
ax = fig.add_subplot(gsp[0, 0])
M = K.loc[SRC, SRC].to_numpy(float)
im = ax.imshow(M, cmap="Blues", vmin=0.3, vmax=1.0, aspect="auto")
for i in range(len(SRC)):
    for j in range(len(SRC)):
        ax.text(j, i, f"{M[i, j]:.2f}".lstrip("0"), ha="center", va="center", fontsize=5.0,
                color="white" if M[i, j] > 0.72 else INK, fontweight="bold" if i == j else "normal")
for b in breaks:
    ax.axhline(b - .5, color="white", lw=1.5); ax.axvline(b - .5, color="white", lw=1.5)
ax.set_xticks(range(len(SRC)), [LAB[s] for s in SRC], rotation=50, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(SRC)), [LAB[s] for s in SRC])
ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02); cb.ax.tick_params(labelsize=5.5, width=.5, length=2)
cb.outline.set_linewidth(.5); cb.set_label("Mean weighted κ", fontsize=6)

# B focal terms (top right, two sub-panels)
sub = gsp[0, 1].subgridspec(1, len(A.focal), wspace=0.62)
axd = []
for j, go in enumerate(A.focal):
    a = fig.add_subplot(sub[0, j]); axd.append(a)
    d = F[F.go == go].set_index("source").reindex(SRC)
    left = np.zeros(len(SRC))
    for k in range(1, 6):
        w = 100 * d[f"n_{k}"] / d.n_runs
        a.barh(range(len(SRC)), w, left=left, color=SEQ[k - 1], edgecolor="white", lw=0.5, height=0.72)
        left += w.to_numpy()
    for i, (_, r) in enumerate(d.iterrows()):
        a.text(103, i, f"{r['mean']:.2f} ({r.rank_in_source})", va="center", ha="left", fontsize=5.2,
               color=INK, clip_on=False)
    a.set_xlim(0, 100); a.set_ylim(len(SRC) - .5, -.5); a.set_xticks([0, 50, 100])
    for b in breaks: a.axhline(b - .5, color=MUTED, lw=.4, ls=(0, (2, 2)))
    for s_ in ("top", "right"): a.spines[s_].set_visible(False)
    t = names.get(go, go); a.set_title(t[0].upper() + t[1:], loc="left", pad=4, fontsize=6.5)
    a.set_xlabel("% of runs")
    if j == 0: a.set_yticks(range(len(SRC)), [LAB[s] for s in SRC])
    else: a.set_yticks(range(len(SRC)), [""] * len(SRC))
axd[0].legend(handles=[Patch(visible=False, label="Score:")] +
                      [Patch(facecolor=SEQ[k - 1], edgecolor=MUTED if k == 1 else SEQ[k - 1], lw=.3, label=str(k))
                       for k in range(1, 6)], ncol=6, frameon=False, handlelength=1, handletextpad=.3,
              columnspacing=.7, loc="upper left", bbox_to_anchor=(-0.05, -0.115), fontsize=5.5)

# C consensus terms with high mean NES: NES in each species and the mean across species
gc = outer[1].subgridspec(1, 3, width_ratios=[0.27, 0.50, 0.23], wspace=0.0)
axc = fig.add_subplot(gc[0, 1])
sp3 = [s for s in SPECIES if s != "Mean"]
y = np.arange(len(TOP))
allm = mnes.to_numpy()
axc.axvspan(np.quantile(allm, .25), np.quantile(allm, .75), color="#f0f0f0", zorder=0, lw=0)
axc.axvline(np.median(allm), color=MUTED, lw=.5, zorder=1)
axc.axvline(cut, color=MUTED, lw=.6, ls=(0, (3, 2)), zorder=1)
for i, r in TOP.iterrows():
    v = [r[f"NES_{s}"] for s in sp3]
    axc.plot([min(v), max(v)], [i, i], color="#b0b0b0", lw=.8, zorder=2, solid_capstyle="butt")
    for s in sp3:
        axc.scatter(r[f"NES_{s}"], i, s=20, marker=SPM[s], facecolor=SPC[s], edgecolor="white", lw=.3, zorder=3)
    axc.scatter(r["NES_Mean"], i, s=34, marker="|", color=INK, lw=1.4, zorder=4)
lo = min(TOP[[f"NES_{s}" for s in sp3]].min().min(), np.quantile(allm, .25)) if len(TOP) else 0.5
hi = TOP[[f"NES_{s}" for s in sp3]].max().max() if len(TOP) else 2
axc.set_xlim(np.floor(lo * 10) / 10 - .05, np.ceil(hi * 10) / 10 + .05)
axc.set_ylim(len(TOP) - .5, -.7)
lab = [abbrev(t) for t in TOP.term_name]
axc.set_yticks(y, lab); axc.tick_params(axis="y", length=0, pad=4, labelsize=7)
axc.set_xlabel("Normalized enrichment score (NES)")
axc.grid(axis="x", color=GRID, lw=.4); axc.set_axisbelow(True)
for s_ in ("top", "right", "left"): axc.spines[s_].set_visible(False)
# reference-line labels
axc.annotate("Median of all terms\n(shaded: middle 50%)", xy=(np.median(allm), 1), xycoords=("data", "axes fraction"),
             xytext=(0, 2), textcoords="offset points", ha="center", va="bottom", fontsize=4.8, color=MUTED)
axc.annotate(f"Top {A.nes_top:.0%} of\nall terms", xy=(cut, 1), xycoords=("data", "axes fraction"),
             xytext=(0, 2), textcoords="offset points", ha="center", va="bottom", fontsize=4.8, color=MUTED)
# right-hand columns: rank by mean NES, mean LLM score, genes
cols = [("Rank\n(mean NES)", lambda r: f"{int(r.rank_mean_NES)}"), ("Mean\nscore", lambda r: f"{r.mean_score_all_runs:.2f}"),
        ("Genes", lambda r: f"{int(r['size'])}" if not pd.isna(r['size']) else "")]
tx = axc.get_yaxis_transform()
for c, (hdr, fn) in enumerate(cols):
    xc = 1.08 + 0.15 * c
    axc.text(xc, -0.95, hdr, transform=tx, ha="center", va="bottom", fontsize=6.5, color=INK, clip_on=False)
    for i, r in TOP.iterrows():
        axc.text(xc, i, fn(r), transform=tx, ha="center", va="center", fontsize=6.5, clip_on=False)
axc.legend(handles=[Line2D([], [], marker=SPM[s], ls="", mfc=SPC[s], mec="white", mew=.3, ms=4.5, label=FULL.get(s, s)) for s in sp3] +
                   [Line2D([], [], marker="|", ls="", color=INK, mew=1.4, ms=6, label="Mean")],
           ncol=len(sp3) + 1, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 0.0), bbox_transform=fig.transFigure,
           fontsize=6.5, handletextpad=.2, columnspacing=1.4)

fig.canvas.draw()
letter_title(fig, 0.012, ax.get_position().y1 + 0.030, "A", "Run-vs-run agreement")
letter_title(fig, 0.505, axd[0].get_position().y1 + 0.030, "B", "Scores for two focal terms")
letter_title(fig, 0.012, axc.get_position().y1 + 0.045, "C", "Consensus terms with the highest enrichment")
base = os.path.join(A.out, "Fig4_LLM_consistency")
save(fig, base); plt.close(fig)

# ------------------------------------------------------------------ caption (Science style) with current counts
focal_txt = " and ".join(names.get(g, g) for g in A.focal)
cap = f"""**Fig. 4. Independent LLM raters agree on which GO terms are relevant to a response to *Philornis downsi*.** \
(**A**) Mean weighted κ between runs. Runs are labeled with the model followed by account type (I, informed account; \
N, naive account; a and b, independent runs on the same account and model). (**B**) Scores given to {focal_txt}; \
numbers show mean score followed by rank among {len(ids):,} terms in parentheses. (**C**) Consensus terms (median score \
≥{A.consensus_score} across all {int(meta.n_runs.sum())} runs; {len(Cdf)} terms) whose NES averaged across species is in \
the top {A.nes_top:.0%} of all {len(ids):,} tested terms ({len(TOP)} terms). Symbols, NES in each species; black ticks, \
mean across species. The solid line and shading show the median and middle 50% of mean NES across all tested terms; the \
dashed line marks the top {A.nes_top:.0%}. Rank, rank by mean NES; mean score, mean LLM score across all runs; genes, \
gene-set size."""
open(base + "_caption.md", "w").write(cap + "\n")
open(base + "_caption.txt", "w").write(cap.replace("**", "").replace("*", "") + "\n")
print("wrote", base, "(.pdf/.png/.tiff, caption) and panel source data in", A.out)
