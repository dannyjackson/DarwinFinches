#!/usr/bin/env python3
"""
Validate every run of one source: each scores.tsv must have exactly one row per GO ID in all_pathways.txt, no extra
IDs, and integer scores from 1 to 5. Prints one line per run and exits with status 1 if any run fails.
A failed run should be re-run from scratch by a fresh rater (STEP 2 of prompts/coordinating_prompt.txt).

Usage:
  export PROJDIR=/path/to/project
  python3 ${PROJDIR}/LLM_evaluation/scripts/validate_runs.py --source informed.Opus5_5.run_a
"""
import argparse, glob, os, sys
import pandas as pd

ROOT = os.path.join(os.environ["PROJDIR"], "LLM_evaluation") if "PROJDIR" in os.environ else \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--source", required=True)
ap.add_argument("--terms", default=os.path.join(ROOT, "runs", "all_pathways.txt"))
ap.add_argument("--runs-dir", default=os.path.join(ROOT, "runs"))
A = ap.parse_args()

ref = [l.split(" ", 1)[0] for l in open(A.terms).read().splitlines()[1:] if l.strip()]
refset = set(ref)
failed = []
run_dirs = sorted(glob.glob(os.path.join(A.runs_dir, A.source, "run_[0-9]*")))
if not run_dirs: sys.exit(f"no run folders in {os.path.join(A.runs_dir, A.source)}")
for d in run_dirs:
    run, f = os.path.basename(d), os.path.join(d, "scores.tsv")
    try:
        s = pd.read_csv(f, sep="\t", dtype=str)
        s.columns = [c.strip().lower() for c in s.columns]
        ids = s.pathway.str.strip()
        sc = pd.to_numeric(s.score, errors="coerce")
        problems = []
        if len(s) != len(ref): problems.append(f"{len(s)} rows")
        if ids.duplicated().any(): problems.append(f"{int(ids.duplicated().sum())} duplicated IDs")
        if refset - set(ids): problems.append(f"{len(refset - set(ids))} missing IDs")
        if set(ids) - refset: problems.append(f"{len(set(ids) - refset)} extra IDs")
        bad = ~(sc.isin([1, 2, 3, 4, 5]))
        if bad.any(): problems.append(f"{int(bad.sum())} scores not integers 1-5")
    except Exception as e:
        problems = [f"unreadable ({e})"]
    print(run, "OK" if not problems else "FAIL: " + "; ".join(problems))
    if problems: failed.append(run)
print(f"{len(run_dirs) - len(failed)} of {len(run_dirs)} runs passed")
sys.exit(1 if failed else 0)
