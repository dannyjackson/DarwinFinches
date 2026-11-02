#!/usr/bin/env python3
"""
Figure 4: consistency of LLM a priori GO-term relevance scores, and the consensus terms with the highest
enrichment.

  A  Mean ICC(2,1) between the runs of every pair of sources (heatmap; diagonal = pairs of runs within a source)
  B  Agreement between runs: ICC(2,1) for every pair of runs, averaged within groups defined by how the two runs
     differ (same source; same account and model, different batch; same model, different account; same account,
     different model; different model and account; different AI provider)
  C  Score distributions for two focal terms (default: blood coagulation, adaptive immune response)
  E, F  Venn diagrams of the terms whose median score across each model's runs is >= 4 (E) or = 5 (F)
  Panels C-F and the consensus terms leave out the models in --drop-models (default: GPT 6.1 Sol).
  D  Consensus terms (median score >= 4 across all runs pooled) whose NES averaged across species is in the
     top 10% of all tested terms
Also writes table_top_terms.tsv (the --top-n terms with the highest mean NES and their LLM scores) and the ICC tables:
table_ICC_within_sources.tsv (ICC(2,1) and ICC(2,k) of each source's k runs) and table_ICC_between_sources.tsv
(ICC(2,1) and ICC(2,k) across the sources' per-term mean scores).

Inputs (defaults are in the repository, ${SCRIPTDIR}; every path can be overridden):
  ${SCRIPTDIR}/data/matrices/score_matrix.<informed|naive>.[Provider.]<model>.<batch|account>.tsv
        e.g. score_matrix.informed.Opus5_5.run_a.tsv, score_matrix.naive.ChatGPT.gpt6_1sol.N1.tsv
        columns: pathway, term_name, run_001 ... run_N  (scores 1-5), one matrix per source
        informed batches run_a / run_b (or run2) are labelled Ia / Ib; naive accounts are labelled with the code in
        the file name (N1, N2). N1 and N2 are the naive accounts for that model, not the same people across models.
  ${SCRIPTDIR}/data/accounts.tsv  which account each source ran on (columns source, account); used to group pairs of
        runs in panel B. Sources on the same account share an account ID.
  ${SCRIPTDIR}/data/{cra,for,par}.fst.gsea.tsv  fgsea tables (columns pathway, NES, padj, size, ...),
        shown as VF (vegetarian finch), MGF (medium ground finch), STF (small tree finch)

Usage:
  export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
  export SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
  cd $PROJDIR
  python3 $SCRIPTDIR/scripts/make_figure4.py

Outputs (in --out, default ${PROJDIR}/LLM_evaluation/results):
  Fig4_LLM_consistency.pdf (vector, fonts embedded), .png and .tiff (600 dpi), caption (.md, .txt)
  panelA_ICC_matrix.tsv, panelB_ICC_groups.tsv, panelB_ICC_source_pairs.tsv, panelC_focal_terms.tsv,
  panelD_consensus_terms_NES.tsv,
  table_venn_terms.tsv (every term in panels E and F: which models scored it >= 4 or 5, and its NES, P, adjusted P
  and NES rank in each species and for the mean), table_ICC_within_sources.tsv, table_ICC_between_sources.tsv, table_top_terms.tsv,
  table_ttest_term_sets_NES.tsv (one-sided Welch's t-tests of NES, each a priori term set vs all other terms, for
  each species and the mean NES), sources.tsv

Requires: Python 3, numpy, pandas, matplotlib (Pillow for TIFF).
"""
import argparse, glob, hashlib, os, re, sys, warnings
import numpy as np, pandas as pd
from scipy.stats import ttest_ind
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
warnings.filterwarnings("ignore")

# Repository (${SCRIPTDIR}): scripts/, prompts/ (rater prompt, term list), data/ (fgsea tables), data/matrices/.
# Results are written outside the repository, to ${PROJDIR}/LLM_evaluation/results/.
CODE = os.environ.get("SCRIPTDIR") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROMPTS, DATA = os.path.join(CODE, "prompts"), os.path.join(CODE, "data")
WORK = os.path.join(os.environ.get("PROJDIR", os.getcwd()), "LLM_evaluation")
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--matrices", default=os.path.join(DATA, "matrices"))
ap.add_argument("--gsea", nargs="+", default=[f"{lab}={os.path.join(DATA, f + '.fst.gsea.tsv')}"
                                               for lab, f in (("VF", "cra"), ("MGF", "for"), ("STF", "par"))],
                help="LABEL=path for each species; the label is what the figure shows")
ap.add_argument("--accounts", default=os.path.join(DATA, "accounts.tsv"),
                help="TSV with columns source and account: which account each source ran on")
ap.add_argument("--out", default=os.path.join(WORK, "results"))
ap.add_argument("--focal", nargs="+", default=["GO:0007596", "GO:0002250"])
ap.add_argument("--names", nargs="*", default=[], help="LABEL=full name overrides for the species legend, e.g. VF='Vegetarian finch'")
ap.add_argument("--consensus-score", type=int, default=4, help="main panel D: score threshold for consensus terms (all runs pooled)")
ap.add_argument("--nes-top", type=float, default=0.10, help="main panel D: show consensus terms whose mean NES is in this top fraction of all terms")
ap.add_argument("--drop-models", nargs="*", default=["gpt6_1sol"],
                help="models left out of the consensus terms, panels C-F and the t-test (default: GPT 6.1 Sol, the model "
                     "with the lowest within-source ICC); panels A and B always show every source")
ap.add_argument("--top-n", type=int, default=10, help="number of top terms by mean NES listed in table_top_terms.tsv")
A = ap.parse_args()
os.makedirs(A.out, exist_ok=True)

# ------------------------------------------------------------------ load matrices
files = sorted(glob.glob(os.path.join(A.matrices, "score_matrix.*.tsv")))
if not files: sys.exit(f"No score_matrix.*.tsv files in {A.matrices}")
# File names: score_matrix.<informed|naive>.[Provider.]<model>.<batch or naive account code>.tsv
#   informed: last token = run_a / run_b / run2 (batch a or b on the first author's accounts)
#   naive:    last token = the naive account code for that model (N1, N2), used as the label
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
        print(f"WARNING: {b} is byte-identical to {seen[h]}; skipping it (check that the right file was copied)")
        continue
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
    meta.append(dict(source=b, model=model.lower(), naive=naive, account=last if naive else "first author",
                     batch="" if naive else batch_letter(last), n_runs=len(runs)))
    mats[b] = X.to_numpy(float)
meta = pd.DataFrame(meta)
meta["code"] = [h if n else f"I{bt}" for h, n, bt in zip(meta.account, meta.naive, meta.batch)]
dup = meta[meta.duplicated(["model", "code"], keep=False)]
if len(dup): sys.exit("two matrices have the same model and account code:\n" + "\n".join(dup.source))
meta["label"] = [f"{pretty_model(m)} · {c}" for m, c in zip(meta.model, meta.code)]
meta["order"] = [(model_rank(m), n, c) for m, n, c in zip(meta.model, meta.naive, meta.code)]
meta = meta.sort_values("order").drop(columns="order").reset_index(drop=True)
PROVIDER = {"opus5_5": "Anthropic", "fable5_1": "Anthropic", "gpt6_1sol": "OpenAI", "gpt6astra": "OpenAI"}
meta["provider"] = [PROVIDER.get(m, "OpenAI" if m.startswith("gpt") else "Anthropic") for m in meta.model]
if not os.path.exists(A.accounts): sys.exit(f"account table not found: {A.accounts}")
acc = pd.read_csv(A.accounts, sep="\t", comment="#", dtype=str)
acc = dict(zip(acc.source.str.strip().str.lower(), acc.account.str.strip()))   # source names matched ignoring case
missing = [x for x in meta.source if x.lower() not in acc]
if missing: sys.exit(f"{A.accounts} has no account for: " + ", ".join(missing))
meta["account_id"] = [acc[x.lower()] for x in meta.source]
SRC = meta.source.tolist(); LAB = dict(zip(meta.source, meta.label))
meta.to_csv(os.path.join(A.out, "sources.tsv"), sep="\t", index=False)
print(f"{len(SRC)} sources, {meta.n_runs.sum()} runs")

# ------------------------------------------------------------------ agreement (ICC)
# ICC(2,*) (Shrout & Fleiss 1979): two-way random effects (terms = targets, runs = random raters), absolute agreement.
def icc2(D):
    """ICC(2,1) and ICC(2,k) for an n x k matrix (rows = targets, columns = raters)."""
    D = np.asarray(D, float); n, k = D.shape; m = D.mean()
    ssr = k * ((D.mean(1) - m) ** 2).sum(); ssc = n * ((D.mean(0) - m) ** 2).sum()
    sse = ((D - m) ** 2).sum() - ssr - ssc
    msr, msc, mse = ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))
    return ((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n),
            (msr - mse) / (msr + (msc - mse) / n))

def icc21_pairs(X, Y):
    """ICC(2,1) for every pair of runs (a column of X with a column of Y): icc2() with k = 2, vectorized."""
    n = X.shape[0]
    Sx, Sy = X.sum(0)[:, None], Y.sum(0)[None]
    Sxx, Syy = (X ** 2).sum(0)[:, None], (Y ** 2).sum(0)[None]
    Sxy = X.T @ Y
    m = (Sx + Sy) / (2 * n)                                   # grand mean of the pair
    sst = Sxx + Syy - 2 * n * m ** 2                           # total sum of squares
    ssr = (Sxx + Syy + 2 * Sxy) / 2 - 2 * n * m ** 2           # between terms
    ssc = n * ((Sx / n - m) ** 2 + (Sy / n - m) ** 2)          # between the two runs
    msr, msc, mse = ssr / (n - 1), ssc, (sst - ssr - ssc) / (n - 1)
    return (msr - mse) / (msr + mse + 2 * (msc - mse) / n)

# within each source: reliability of one run, ICC(2,1), and of the mean of its k runs, ICC(2,k)
W_ = pd.DataFrame([dict(source=x, label=LAB[x], n_runs=mats[x].shape[1],
                        **dict(zip(["ICC_2_1", "ICC_2_k"], icc2(mats[x])))) for x in SRC])
W_.round(4).to_csv(os.path.join(A.out, "table_ICC_within_sources.tsv"), sep="\t", index=False)
# between sources: the sources' per-term mean scores as raters
MEANS = np.column_stack([mats[x].mean(1) for x in SRC])
b1, bk = icc2(MEANS)
pd.DataFrame([dict(n_sources=len(SRC), n_terms=len(ids), ICC_2_1=b1, ICC_2_k=bk)]).round(4).to_csv(
    os.path.join(A.out, "table_ICC_between_sources.tsv"), sep="\t", index=False)
print(f"ICC(2,1) between the {len(SRC)} sources' mean scores: {b1:.3f}")

# every pair of runs, grouped by how the two runs differ
GROUPS = ["Same source", "Same account and model, different batch", "Same model, different account",
          "Same account, different model", "Different model and account", "Different AI provider"]
def pair_group(a, b):
    ma, mb = meta.set_index("source").loc[a], meta.set_index("source").loc[b]
    if a == b: return GROUPS[0]
    if ma.provider != mb.provider: return GROUPS[5]
    same_acc = ma.account_id == mb.account_id
    if ma.model == mb.model: return GROUPS[1] if same_acc else GROUPS[2]
    return GROUPS[3] if same_acc else GROUPS[4]
SP = []
for i_, a in enumerate(SRC):
    for b in SRC[i_:]:
        P = icc21_pairs(mats[a], mats[b])
        v = P[np.triu_indices(P.shape[0], 1)] if a == b else P.ravel()
        SP.append(dict(source_1=a, source_2=b, label_1=LAB[a], label_2=LAB[b], group=pair_group(a, b),
                       n_run_pairs=v.size, mean_ICC=v.mean(), sum_ICC=v.sum()))
SP = pd.DataFrame(SP)
SP.drop(columns="sum_ICC").round(4).to_csv(os.path.join(A.out, "panelB_ICC_source_pairs.tsv"), sep="\t", index=False)
GR = (SP.groupby("group").agg(n_source_pairs=("mean_ICC", "size"), n_run_pairs=("n_run_pairs", "sum"),
                              sum_ICC=("sum_ICC", "sum"), min_source_pair=("mean_ICC", "min"),
                              max_source_pair=("mean_ICC", "max"))
      .reindex([g for g in GROUPS if g in set(SP.group)]))
GR.insert(2, "mean_ICC", GR.sum_ICC / GR.n_run_pairs)              # mean over all pairs of runs in the group
GR = GR.drop(columns="sum_ICC").reset_index()
GR.round(4).to_csv(os.path.join(A.out, "panelB_ICC_groups.tsv"), sep="\t", index=False)
print(GR.round(3).to_string(index=False))

# ------------------------------------------------------------------ NES from the fgsea tables
gs = dict(g.split("=", 1) for g in A.gsea)
for sp, f in gs.items():
    if not os.path.exists(f): sys.exit(f"GSEA table for {sp} not found: {f} (set --gsea SP=path ...)")
G = {sp: pd.read_csv(f, sep="\t").set_index("pathway").reindex(ids) for sp, f in gs.items()}
NES = pd.DataFrame({sp: g.NES for sp, g in G.items()})
if NES.isna().any().any(): sys.exit("Some GO IDs are missing from a GSEA table")
NES["Mean"] = NES.mean(axis=1)
SPECIES = list(NES.columns)

DROP = [m.lower() for m in A.drop_models]
KEEP = [x for x in SRC if meta.set_index("source").model[x] not in DROP]   # sources used from here on
if not KEEP: sys.exit("--drop-models leaves no sources")
KEEP_MODELS = [m for m in MODEL_ORDER + sorted(set(meta.model) - set(MODEL_ORDER)) if m in set(meta.model) and m not in DROP]
print(f"consensus terms and panels C-F use {len(KEEP)} sources ({', '.join(pretty_model(m) for m in KEEP_MODELS)}); "
      f"dropped: {', '.join(pretty_model(m) for m in DROP if m in set(meta.model)) or 'none'}")
ALLX = np.hstack([mats[s] for s in KEEP])   # runs of the kept models pooled: terms x runs

# ------------------------------------------------------------------ panel C: focal terms
F = []
for s in KEEP:
    X = mats[s]; mean = X.mean(1); rank = pd.Series(mean).rank(ascending=False, method="min").to_numpy()
    for go in A.focal:
        i = ids.index(go); x = X[i]
        F.append(dict(source=s, label=LAB[s], go=go, term=names.get(go, go), n_runs=len(x),
                      **{f"n_{k}": int((x == k).sum()) for k in range(1, 6)},
                      mean=x.mean(), rank_in_source=int(rank[i])))
F = pd.DataFrame(F); F.round(3).to_csv(os.path.join(A.out, "panelC_focal_terms.tsv"), sep="\t", index=False)


# ------------------------------------------------------------------ main panel D: consensus terms with high NES
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
Cdf.round(4).to_csv(os.path.join(A.out, "panelD_consensus_terms_NES.tsv"), sep="\t", index=False)
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
print(f"panel D: {len(Cdf)} consensus terms, {len(TOP)} with mean NES in the top {A.nes_top:.0%} (>= {cut:.2f})")

# ------------------------------------------------------------------ panels E, F: overlap of the models' high-scoring terms
# For each kept model, the terms whose median score across all of that model's runs is >= 4 (E) or = 5 (F).
MODEL_X = {m: np.hstack([mats[x] for x in KEEP if meta.set_index("source").model[x] == m]) for m in KEEP_MODELS}
VSETS = {t: {m: set(np.array(ids)[np.median(MODEL_X[m], 1) >= t]) for m in KEEP_MODELS} for t in (4, 5)}
# one row per term in either Venn diagram: which models scored it >= 4 (and 5), and its enrichment in each species
VT = sorted(set().union(*VSETS[4].values(), *VSETS[5].values()), key=lambda g: ids.index(g))
rank_mean = NES["Mean"].rank(ascending=False, method="min").astype(int)
rows = []
for go in VT:
    ge4 = [m for m in KEEP_MODELS if go in VSETS[4][m]]; eq5 = [m for m in KEEP_MODELS if go in VSETS[5][m]]
    row = dict(pathway=go, term_name=names[go], size=int(G[SPECIES[0]].loc[go, "size"]) if "size" in G[SPECIES[0]] else np.nan,
               models_ge4=", ".join(pretty_model(m) for m in ge4) or "none", n_models_ge4=len(ge4),
               all_models_ge4=len(ge4) == len(KEEP_MODELS),
               models_5=", ".join(pretty_model(m) for m in eq5) or "none", n_models_5=len(eq5),
               **{f"median_score_{short_model(m)}": float(np.median(MODEL_X[m][ids.index(go)])) for m in KEEP_MODELS})
    for sp in SPECIES[:-1]:
        g_ = G[sp].loc[go]
        row.update({f"NES_{sp}": g_.NES, f"pval_{sp}": g_.get("pval", np.nan), f"padj_{sp}": g_.get("padj", np.nan),
                    f"rank_NES_{sp}": int(NES[sp].rank(ascending=False, method="min")[go])})
    row.update(NES_mean=NES.loc[go, "Mean"], rank_NES_mean=int(rank_mean[go]))
    rows.append(row)
VTAB = pd.DataFrame(rows).sort_values(["n_models_ge4", "rank_NES_mean"], ascending=[False, True])
VTAB.round(4).to_csv(os.path.join(A.out, "table_venn_terms.tsv"), sep="\t", index=False)
for t, sets in VSETS.items():
    print(f"score {'>=4' if t == 4 else '=5'}: " + ", ".join(f"{pretty_model(m)} {len(v)}" for m, v in sets.items())
          + f"; all {len(KEEP_MODELS)} models {len(set.intersection(*sets.values()))}")

# ------------------------------------------------------------------ t-tests: NES of a priori term sets vs all other terms
# One-sided Welch's t-test (the set has higher NES than the other tested terms), for each species and the mean NES.
#   shared_ge4: median score >= 4 within each kept model's runs, in every kept model (centre of panel E)
#   pooled_ge4: median score >= 4 across all runs of the kept models pooled (the consensus terms of panel D)
TSETS = [("shared_ge4", f"Median >=4 in each of the {len(KEEP_MODELS)} models", np.isin(ids, list(set.intersection(*VSETS[4].values())))),
         ("pooled_ge4", "Median >=4 across all runs pooled", np.median(ALLX, 1) >= 4)]
TT = []
for key, desc, m_ in TSETS:
    for sp in SPECIES:
        v = NES[sp].to_numpy()
        r = ttest_ind(v[m_], v[~m_], equal_var=False, alternative="greater") if m_.sum() >= 2 else None
        TT.append(dict(set=key, definition=desc, species=sp, n_set=int(m_.sum()), n_other=int((~m_).sum()),
                       mean_NES_set=v[m_].mean() if m_.any() else np.nan, mean_NES_other=v[~m_].mean(),
                       difference=(v[m_].mean() - v[~m_].mean()) if m_.any() else np.nan,
                       t=r.statistic if r else np.nan, df=r.df if r else np.nan, p_one_sided=r.pvalue if r else np.nan,
                       terms="; ".join(names[g] for g in np.array(ids)[m_]) if sp == SPECIES[0] else ""))
TT = pd.DataFrame(TT)
TT.round(4).to_csv(os.path.join(A.out, "table_ttest_term_sets_NES.tsv"), sep="\t", index=False)
print(TT.drop(columns=["terms", "definition"]).round(3).to_string(index=False))

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

# ================================================================== MAIN FIGURE: A, B (top row), C, D (bottom row)
fig = plt.figure(figsize=(W, W * 1.34))
outer = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.66, 0.40], hspace=0.50, left=0.105, right=0.945, top=0.955, bottom=0.015)
gsp = outer[0].subgridspec(1, 4, width_ratios=[0.02, 1.0, 0.44, 0.58], wspace=0.0)   # [| A | gap + B labels | B]
bot = outer[1].subgridspec(1, 2, width_ratios=[0.40, 0.60], wspace=0.0)               # [C | D]
model_of = meta.set_index("source").model.to_dict()
breaks = [i for i in range(1, len(SRC)) if model_of[SRC[i]] != model_of[SRC[i - 1]]]
breaksK = [i for i in range(1, len(KEEP)) if model_of[KEEP[i]] != model_of[KEEP[i - 1]]]

# B agreement by how two runs differ: grey points = one per pair of sources (upper half of each row);
#   black diamond = mean over all pairs of runs in the group (offset below the points so it does not hide them)
ax = fig.add_subplot(gsp[0, 3])
ylab = {GROUPS[0]: "Same source", GROUPS[1]: "Same account and model,\ndifferent batch",
        GROUPS[2]: "Same model,\ndifferent account", GROUPS[3]: "Same account,\ndifferent model",
        GROUPS[4]: "Different model\nand account", GROUPS[5]: "Different AI provider"}
rng = np.random.default_rng(1)
for i, r in GR.iterrows():
    v = SP[SP.group == r.group].mean_ICC.to_numpy()
    ax.scatter(v, i - 0.14 + rng.uniform(-0.11, 0.11, v.size), s=7, color="#9a9a9a", lw=0, alpha=.85, zorder=2)
    ax.scatter(r.mean_ICC, i + 0.22, s=26, marker="D", color=INK, zorder=3)
    ax.text(1.03, i + 0.22, f"{r.mean_ICC:.2f}", transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=6)
ax.set_yticks(range(len(GR)), [ylab[g] for g in GR.group]); ax.set_ylim(len(GR) - .45, -.65)
ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1], ["0", "0.5", "1"]); ax.set_xticks([.25, .75], minor=True)
ax.set_xlabel("ICC(2,1) between\npairs of runs")
ax.grid(axis="x", which="minor", color=GRID, lw=.4)
ax.grid(axis="x", color=GRID, lw=.4); ax.set_axisbelow(True); ax.tick_params(axis="y", length=0)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.text(1.03, -0.62, "Mean", transform=ax.get_yaxis_transform(), va="bottom", ha="left", fontsize=6)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#9a9a9a", ms=3, label="Pair of sources"),
                   Line2D([], [], marker="D", ls="", color=INK, ms=3.6, label="All pairs of runs")],
          frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=1, fontsize=6, handletextpad=.2)

# A heatmap: mean ICC(2,1) between the runs of each pair of sources (diagonal = pairs of runs within a source)
axh = fig.add_subplot(gsp[0, 1])
K = pd.DataFrame(np.nan, index=SRC, columns=SRC)
for r in SP.itertuples():
    K.loc[r.source_1, r.source_2] = K.loc[r.source_2, r.source_1] = r.mean_ICC
K.round(4).to_csv(os.path.join(A.out, "panelA_ICC_matrix.tsv"), sep="\t")
M = K.to_numpy(float)
im = axh.imshow(M, cmap="Blues", vmin=0.3, vmax=1.0, aspect="auto")
for i in range(len(SRC)):
    for j in range(len(SRC)):
        axh.text(j, i, f"{M[i, j]:.2f}".lstrip("0"), ha="center", va="center", fontsize=4.3,
                 color="white" if M[i, j] > 0.72 else INK, fontweight="bold" if i == j else "normal")
for b in breaks:
    axh.axhline(b - .5, color="white", lw=1.2); axh.axvline(b - .5, color="white", lw=1.2)
axh.set_xticks(range(len(SRC)), [LAB[x] for x in SRC], rotation=55, ha="right", rotation_mode="anchor", fontsize=5.3)
axh.set_yticks(range(len(SRC)), [LAB[x] for x in SRC], fontsize=5.3)
axh.tick_params(length=0)
for s_ in axh.spines.values(): s_.set_visible(False)
cb = fig.colorbar(im, ax=axh, fraction=0.045, pad=0.02); cb.ax.tick_params(labelsize=5.5, width=.5, length=2)
cb.outline.set_linewidth(.5); cb.set_label("Mean ICC(2,1)", fontsize=6)

# C focal terms (bottom left, two narrow sub-panels)
sub = bot[0, 0].subgridspec(1, 5, width_ratios=[0.40, 0.55, 0.36, 0.55, 0.30], wspace=0.0)
axd = []
for j, go in enumerate(A.focal):
    a = fig.add_subplot(sub[0, 1 + 2 * j]); axd.append(a)
    d = F[F.go == go].set_index("source").reindex(KEEP)
    left = np.zeros(len(KEEP))
    for k in range(1, 6):
        w = 100 * d[f"n_{k}"] / d.n_runs
        a.barh(range(len(KEEP)), w, left=left, color=SEQ[k - 1], edgecolor="white", lw=0.5, height=0.72)
        left += w.to_numpy()
    for i, (_, r) in enumerate(d.iterrows()):
        a.text(106, i, f"{r['mean']:.2f} ({r.rank_in_source})", va="center", ha="left", fontsize=4.8,
               color=INK, clip_on=False)
    a.set_xlim(0, 100); a.set_ylim(len(KEEP) - .5, -.5); a.set_xticks([0, 50, 100])
    for b in breaksK: a.axhline(b - .5, color=MUTED, lw=.4, ls=(0, (2, 2)))
    for s_ in ("top", "right"): a.spines[s_].set_visible(False)
    t = names.get(go, go); t = t[0].upper() + t[1:]
    a.set_title(t.replace(" response", "\nresponse") if len(t) > 20 else t, loc="left", pad=3, fontsize=6)
    a.set_xlabel("% of runs"); a.tick_params(axis="x", labelsize=5.5)
    if j == 0: a.set_yticks(range(len(KEEP)), [LAB[s] for s in KEEP], fontsize=5.3)
    else: a.set_yticks(range(len(KEEP)), [""] * len(KEEP))
axd[0].legend(handles=[Patch(visible=False, label="Score:")] +
                      [Patch(facecolor=SEQ[k - 1], edgecolor=MUTED if k == 1 else SEQ[k - 1], lw=.3, label=str(k))
                       for k in range(1, 6)], ncol=6, frameon=False, handlelength=1, handletextpad=.3,
              columnspacing=.6, loc="upper left", bbox_to_anchor=(-0.75, -0.15), fontsize=5.5)

# D consensus terms with high mean NES: NES in each species and the mean across species
gc = bot[0, 1].subgridspec(3, 3, width_ratios=[0.37, 0.37, 0.26], height_ratios=[0.10, 1, 0.20], wspace=0.0, hspace=0.0)
axc = fig.add_subplot(gc[1, 1])
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
import textwrap
lab = ["\n".join(textwrap.wrap(abbrev(t), 26)) for t in TOP.term_name]
axc.set_yticks(y, lab); axc.tick_params(axis="y", length=0, pad=4, labelsize=6.3)
axc.set_xlabel("NES"); axc.tick_params(axis="x", labelsize=5.5)
axc.grid(axis="x", color=GRID, lw=.4); axc.set_axisbelow(True)
for s_ in ("top", "right", "left"): axc.spines[s_].set_visible(False)
# reference-line labels
axc.annotate("Median", xy=(np.median(allm), 1), xycoords=("data", "axes fraction"),
             xytext=(-1, 2), textcoords="offset points", ha="right", va="bottom", fontsize=5, color=MUTED)
axc.annotate(f"Top {A.nes_top:.0%}", xy=(cut, 1), xycoords=("data", "axes fraction"),
             xytext=(1, 2), textcoords="offset points", ha="left", va="bottom", fontsize=5, color=MUTED)
# right-hand columns: rank by mean NES, mean LLM score, genes
cols = [("NES\nrank", lambda r: f"{int(r.rank_mean_NES)}"), ("Mean\nscore", lambda r: f"{r.mean_score_all_runs:.2f}"),
        ("Genes", lambda r: f"{int(r['size'])}" if not pd.isna(r['size']) else "")]
tx = axc.get_yaxis_transform()
for c, (hdr, fn) in enumerate(cols):
    xc = 1.16 + 0.24 * c
    axc.text(xc, -0.95, hdr, transform=tx, ha="center", va="bottom", fontsize=6, color=INK, clip_on=False)
    for i, r in TOP.iterrows():
        axc.text(xc, i, fn(r), transform=tx, ha="center", va="center", fontsize=6, clip_on=False)
axc.legend(handles=[Line2D([], [], marker=SPM[s], ls="", mfc=SPC[s], mec="white", mew=.3, ms=4.5, label=FULL.get(s, s)) for s in sp3] +
                   [Line2D([], [], marker="|", ls="", color=INK, mew=1.4, ms=6, label="Mean")],
           ncol=2, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), fontsize=6, handletextpad=.2,
           columnspacing=1.0)

# E, F Venn diagrams: terms whose median score across a model's runs is >= 4 (E) or = 5 (F)
VCOL = {"opus5_5": "#d0743c", "fable5_1": "#8c6bb1", "gpt6_1sol": "#7fb800", "gpt6astra": "#2b8cbe"}
def venn(ax, sets):
    ms = list(sets)
    ax.set_aspect("equal"); ax.axis("off")
    if len(ms) == 3:
        C = [(-0.62, 0.36), (0.62, 0.36), (0.0, -0.70)]
        lab_xy = [(-1.45, 1.32, "right"), (1.45, 1.32, "left"), (0.0, -2.08, "center")]
        a, b, c = (sets[m] for m in ms)
        reg = {(-1.05, 0.62): a - b - c, (1.05, 0.62): b - a - c, (0.0, -1.20): c - a - b,
               (0.0, 0.80): (a & b) - c, (-0.70, -0.38): (a & c) - b, (0.70, -0.38): (b & c) - a, (0.0, 0.02): a & b & c}
        ax.set_xlim(-2.3, 2.3); ax.set_ylim(-2.35, 1.65)
    elif len(ms) == 2:
        C = [(-0.55, 0.0), (0.55, 0.0)]; lab_xy = [(-1.3, 1.1, "right"), (1.3, 1.1, "left")]
        a, b = (sets[m] for m in ms)
        reg = {(-1.0, 0.0): a - b, (1.0, 0.0): b - a, (0.0, 0.0): a & b}
        ax.set_xlim(-2.2, 2.2); ax.set_ylim(-1.3, 1.4)
    else:
        ax.text(0.5, 0.5, "Venn diagram needs 2 or 3 models", transform=ax.transAxes, ha="center"); return
    for (x, y), m in zip(C, ms):
        ax.add_patch(plt.Circle((x, y), 1.0, facecolor=VCOL.get(m, "#999999"), alpha=.22, lw=0))
        ax.add_patch(plt.Circle((x, y), 1.0, fill=False, edgecolor=VCOL.get(m, "#999999"), lw=1.0))
    for (x, y, ha), m in zip(lab_xy, ms):
        ax.text(x, y, f"{pretty_model(m)}\n({len(sets[m])})", ha=ha, va="center", fontsize=6.5, color=INK)
    for (x, y), v in reg.items():
        ax.text(x, y, str(len(v)), ha="center", va="center", fontsize=7,
                fontweight="bold" if (x, y) in ((0.0, 0.02), (0.0, 0.0)) else "normal")
gv = outer[2].subgridspec(1, 4, width_ratios=[0.08, 1, 0.12, 1], wspace=0.0)
axv = [fig.add_subplot(gv[0, 1]), fig.add_subplot(gv[0, 3])]
venn(axv[0], VSETS[4]); venn(axv[1], VSETS[5])

fig.canvas.draw()
letter_title(fig, 0.012, axh.get_position().y1 + 0.022, "A", "Agreement between sources")
letter_title(fig, 0.575, axh.get_position().y1 + 0.022, "B", "Agreement by how runs differ")
letter_title(fig, 0.012, axd[0].get_position().y1 + 0.045, "C", "Scores for two focal terms")
letter_title(fig, 0.44, axd[0].get_position().y1 + 0.045, "D", "Consensus terms with the highest enrichment")
letter_title(fig, 0.012, axv[0].get_position().y1 + 0.004, "E", "Terms scored ≥4 by each model")
letter_title(fig, 0.51, axv[1].get_position().y1 + 0.004, "F", "Terms scored 5 by each model")
base = os.path.join(A.out, "Fig4_LLM_consistency")
save(fig, base); plt.close(fig)

# ------------------------------------------------------------------ caption (Science style) with current counts
focal_txt = " and ".join(names.get(g, g) for g in A.focal)
kept_names = [pretty_model(m) for m in KEEP_MODELS]
models_txt = ", ".join(kept_names[:-1]) + " and " + kept_names[-1] if len(kept_names) > 1 else kept_names[0]
dropped = [pretty_model(m) for m in DROP if m in set(meta.model)]
drop_txt = (f"Panels C to F exclude {' and '.join(dropped)}, the model with the lowest within-source agreement. "
            if dropped else "")
cap = f"""**Fig. 4. Independent LLM raters agree on which GO terms are relevant to a response to *Philornis downsi*.** \
(**A**) Mean ICC(2,1) between the runs of each pair of sources; the diagonal is the mean over pairs of runs within a \
source. Sources are labeled with the model followed by account type (I, informed account; N1 and N2, naive accounts, \
numbered separately for each model; a and b, independent batches on the same account and model). (**B**) Agreement \
between pairs of runs (ICC(2,1)), grouped by how the two runs differ. Grey points, mean for each pair of sources; black \
diamonds and numbers, mean over all pairs of runs in the group. {drop_txt}(**C**) Scores given to {focal_txt} by each \
source; numbers show mean score followed by rank among {len(ids):,} terms in parentheses. (**D**) Consensus terms (median \
score ≥{A.consensus_score} across all {ALLX.shape[1]} runs of {models_txt}; {len(Cdf)} terms) whose NES averaged across species is in \
the top {A.nes_top:.0%} of all {len(ids):,} tested terms ({len(TOP)} terms). Symbols, NES in each species; black ticks, \
mean across species. The solid line and shading show the median and middle 50% of mean NES across all tested terms; the \
dashed line marks the top {A.nes_top:.0%}. NES rank, rank by mean NES; mean score, mean LLM score across those runs; genes, \
gene-set size. (**E** and **F**) Number of GO terms whose median score across all runs of each model is ≥4 (E) or 5 (F); \
numbers in parentheses are the totals for each model, and the central number is the terms shared by all {len(KEEP_MODELS)} models."""
open(base + "_caption.md", "w").write(cap + "\n")
open(base + "_caption.txt", "w").write(cap.replace("**", "").replace("*", "") + "\n")
print("wrote", base, "(.pdf/.png/.tiff, caption) and panel source data in", A.out)
