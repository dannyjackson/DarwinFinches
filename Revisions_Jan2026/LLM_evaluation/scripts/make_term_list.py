#!/usr/bin/env python3
"""
Write the list of GO terms given to the LLM raters (all_pathways.txt).

The list is every GO term in the fgsea table, in table order: a header line "pathway term_name", then one line per
term with the GO ID, a space, and the term name. Terms are checked against the other species' tables so that every
listed term was tested in all three species.

Usage:
  export PROJDIR=/path/to/project
  python3 ${PROJDIR}/LLM_evaluation/scripts/make_term_list.py
"""
import argparse, os, sys
import pandas as pd

ROOT = os.path.join(os.environ["PROJDIR"], "LLM_evaluation") if "PROJDIR" in os.environ else \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--gsea", nargs="+", default=[os.path.join(ROOT, "data", "gsea", f + ".fst.gsea.tsv")
                                               for f in ("cra", "for", "par")],
                help="fgsea tables; the first sets the term order")
ap.add_argument("--out", default=os.path.join(ROOT, "runs", "all_pathways.txt"))
A = ap.parse_args()

G = [pd.read_csv(f, sep="\t") for f in A.gsea]
ref = G[0][["pathway", "term_name"]]
for f, g in zip(A.gsea[1:], G[1:]):
    if set(g.pathway) != set(ref.pathway):
        sys.exit(f"{f}: GO terms differ from {A.gsea[0]}")
if ref.pathway.duplicated().any(): sys.exit("duplicate GO IDs in the fgsea table")
os.makedirs(os.path.dirname(A.out), exist_ok=True)
with open(A.out, "w") as fh:
    fh.write("pathway term_name\n")
    for r in ref.itertuples():
        fh.write(f"{r.pathway} {r.term_name}\n")
print(f"wrote {len(ref)} terms to {A.out}")
