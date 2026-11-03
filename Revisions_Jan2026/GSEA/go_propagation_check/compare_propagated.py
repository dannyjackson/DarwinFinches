#!/usr/bin/env python3
"""
Compare the direct-annotation fgsea results (the manuscript) with the propagated-annotation
results: for the focal terms, report term size, NES, nominal P, adjusted P and NES rank in each
species, and rank by mean NES among terms tested in all three species. Also reports overall
rank concordance and how the top 10 by mean NES changed.

Usage:
  python3 compare_propagated.py --direct gsea --propagated gsea_propagated \
      --species cra:VF for:MGF par:STF --focal GO:0002250 GO:0007596 --out gsea_propagated/comparison

Expects <dir>/<species>.fst.gsea.tsv as written by D1_run_fgsea.R (columns pathway, NES, pval,
padj, size, term_name).
"""
import argparse, os
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--direct", required=True)
ap.add_argument("--propagated", required=True)
ap.add_argument("--species", nargs="+", default=["cra:VF", "for:MGF", "par:STF"], help="fileprefix:label")
ap.add_argument("--focal", nargs="+", default=["GO:0002250", "GO:0007596"])
ap.add_argument("--out", required=True, help="output prefix")
A = ap.parse_args()
os.makedirs(os.path.dirname(A.out) or ".", exist_ok=True)
SP = [s.split(":") for s in A.species]

def load(d):
    tabs = {}
    for pre, lab in SP:
        t = pd.read_csv(os.path.join(d, f"{pre}.fst.gsea.tsv"), sep="\t").set_index("pathway")
        t = t[~t.NES.isna()]
        t["rank_NES"] = t.NES.rank(ascending=False, method="min").astype(int)
        tabs[lab] = t
    common = set.intersection(*(set(t.index) for t in tabs.values()))
    mean = pd.DataFrame({lab: t.NES for lab, t in tabs.items()}).loc[sorted(common)]
    mean["mean_NES"] = mean.mean(axis=1)
    mean["rank_mean_NES"] = mean.mean_NES.rank(ascending=False, method="min").astype(int)
    mean["n_positive"] = (mean[[l for _, l in SP]] > 0).sum(axis=1)
    names = pd.concat([t.term_name for t in tabs.values()]).groupby(level=0).first()
    return tabs, mean, names

D, Dm, names = load(A.direct)
P, Pm, namesP = load(A.propagated)
names = names.combine_first(namesP)

rows = []
for go in A.focal:
    for kind, tabs, mean in (("direct", D, Dm), ("propagated", P, Pm)):
        r = dict(term=go, term_name=names.get(go, go), annotation=kind)
        for _, lab in SP:
            t = tabs[lab]
            if go in t.index:
                r.update({f"size_{lab}": int(t.loc[go, "size"]), f"NES_{lab}": round(t.loc[go, "NES"], 3),
                          f"p_{lab}": t.loc[go, "pval"], f"padj_{lab}": round(t.loc[go, "padj"], 3),
                          f"rank_{lab}": int(t.loc[go, "rank_NES"]), f"n_terms_{lab}": len(t)})
            else:
                r.update({f"size_{lab}": np.nan, f"NES_{lab}": np.nan, f"rank_{lab}": np.nan})
        if go in mean.index:
            r.update(mean_NES=round(mean.loc[go, "mean_NES"], 3), rank_mean_NES=int(mean.loc[go, "rank_mean_NES"]),
                     n_terms_common=len(mean))
        rows.append(r)
F = pd.DataFrame(rows)
F.to_csv(A.out + ".focal_terms.tsv", sep="\t", index=False)
print("\n== Focal terms ==")
print(F.to_string(index=False))

# term universe and concordance
common = sorted(set(Dm.index) & set(Pm.index))
rho = spearmanr(Dm.loc[common, "mean_NES"], Pm.loc[common, "mean_NES"]).statistic
print(f"\nterms tested in all three species: direct {len(Dm)}, propagated {len(Pm)}, shared {len(common)}")
print(f"Spearman rho of mean NES across the {len(common)} shared terms: {rho:.3f}")

top = pd.DataFrame({
    "direct_top10": [f"{names.get(g, g)} ({g})" for g in Dm.sort_values('rank_mean_NES').index[:10]],
    "propagated_top10": [f"{names.get(g, g)} ({g})" for g in Pm.sort_values('rank_mean_NES').index[:10]]})
top.to_csv(A.out + ".top10.tsv", sep="\t", index=False)
print("\n== Top 10 by mean NES ==")
print(top.to_string(index=False))

# where the direct top 10 land under propagation
moved = pd.DataFrame({
    "term": [names.get(g, g) for g in Dm.sort_values('rank_mean_NES').index[:10]],
    "rank_direct": Dm.sort_values('rank_mean_NES').rank_mean_NES.values[:10],
    "rank_propagated": [int(Pm.loc[g, "rank_mean_NES"]) if g in Pm.index else -1
                        for g in Dm.sort_values('rank_mean_NES').index[:10]],
    "size_direct": [int(D[SP[0][1]].loc[g, "size"]) for g in Dm.sort_values('rank_mean_NES').index[:10]],
    "size_propagated": [int(P[SP[0][1]].loc[g, "size"]) if g in P[SP[0][1]].index else -1
                        for g in Dm.sort_values('rank_mean_NES').index[:10]]})
moved.to_csv(A.out + ".direct_top10_under_propagation.tsv", sep="\t", index=False)
print("\n== Direct top 10: rank and size under propagation (-1 = not tested, size outside 10-500) ==")
print(moved.to_string(index=False))

# significance
for kind, tabs in (("direct", D), ("propagated", P)):
    sig = {lab: int((t.padj < 0.05).sum()) for lab, t in tabs.items()}
    print(f"\n{kind}: terms with padj < 0.05 per species: {sig}")
print("\nwrote", A.out + ".*.tsv")
