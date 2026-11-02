#!/usr/bin/env python3
"""
Write, or check, the list of GO terms given to the LLM raters (prompts/all_pathways.txt).

The list is every GO term in the fgsea table, in table order: a header line "pathway term_name", then one line per
term with the GO ID, a space, and the term name. Terms are checked against the other species' tables so that every
listed term was tested in all three species. If prompts/all_pathways.txt already exists, the script checks that it
matches and leaves it unchanged.

Usage:
  export PROJDIR=/xdisk/mcnew/finches/dannyjackson/simulations/gsea_power_analysis
  export SCRIPTDIR=${PROJDIR}/DarwinFinches/Revisions_Jan2026/LLM_evaluation/
  cd $PROJDIR
  python3 $SCRIPTDIR/scripts/make_term_list.py
"""
import argparse, os, sys
import pandas as pd

# Repository (${SCRIPTDIR}): scripts/, prompts/ (rater prompt, term list), data/ (fgsea tables), data/matrices/.
CODE = os.environ.get("SCRIPTDIR") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROMPTS, DATA = os.path.join(CODE, "prompts"), os.path.join(CODE, "data")
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--gsea", nargs="+", default=[os.path.join(DATA, f + ".fst.gsea.tsv")
                                               for f in ("cra", "for", "par")],
                help="fgsea tables; the first sets the term order")
ap.add_argument("--out", default=os.path.join(PROMPTS, "all_pathways.txt"))
ap.add_argument("--force", action="store_true", help="overwrite --out if it exists and differs")
A = ap.parse_args()

G = [pd.read_csv(f, sep="\t") for f in A.gsea]
ref = G[0][["pathway", "term_name"]]
for f, g in zip(A.gsea[1:], G[1:]):
    if set(g.pathway) != set(ref.pathway):
        sys.exit(f"{f}: GO terms differ from {A.gsea[0]}")
if ref.pathway.duplicated().any(): sys.exit("duplicate GO IDs in the fgsea table")
text = "pathway term_name\n" + "".join(f"{r.pathway} {r.term_name}\n" for r in ref.itertuples())
if os.path.exists(A.out):
    if open(A.out).read() == text:
        print(f"{A.out} already matches the fgsea tables ({len(ref)} terms); nothing written"); sys.exit(0)
    if not A.force:
        sys.exit(f"{A.out} differs from the list built from the fgsea tables; use --force to overwrite")
os.makedirs(os.path.dirname(A.out), exist_ok=True)
open(A.out, "w").write(text)
print(f"wrote {len(ref)} terms to {A.out}")
