#!/usr/bin/env python3
"""
Compile the validated runs of one source into a score matrix for make_figure4.py (STEP 3 of
prompts/coordinating_prompt.txt).

Writes data/score_matrices/score_matrix.<source>.tsv: one row per GO ID in all_pathways.txt order, with columns
pathway, term_name, run_001 ... run_N; and runs/<source>/run_log.tsv with each run's shuffle seed and the count of
each score 1-5.

Name the source score_matrix.<informed|naive>.[Provider.]<model>.<batch|account holder>, e.g.
informed.Opus5_5.run_a or naive.ChatGPT.gpt6astra.Rebecca, so that make_figure4.py can read the model and account.

Usage:
  export PROJDIR=/path/to/project
  python3 ${PROJDIR}/LLM_evaluation/scripts/validate_runs.py --source informed.Opus5_5.run_a && \
  python3 ${PROJDIR}/LLM_evaluation/scripts/compile_matrix.py --source informed.Opus5_5.run_a
"""
import argparse, glob, os
import pandas as pd

ROOT = os.path.join(os.environ["PROJDIR"], "LLM_evaluation") if "PROJDIR" in os.environ else \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--source", required=True)
ap.add_argument("--terms", default=os.path.join(ROOT, "runs", "all_pathways.txt"))
ap.add_argument("--runs-dir", default=os.path.join(ROOT, "runs"))
ap.add_argument("--out-dir", default=os.path.join(ROOT, "data", "score_matrices"))
A = ap.parse_args()

rows = [l.split(" ", 1) for l in open(A.terms).read().splitlines()[1:] if l.strip()]
M = pd.DataFrame(rows, columns=["pathway", "term_name"]).set_index("pathway")
src_dir = os.path.join(A.runs_dir, A.source)
man = pd.read_csv(os.path.join(src_dir, "run_manifest.tsv"), sep="\t").set_index("run")
log = []
for d in sorted(glob.glob(os.path.join(src_dir, "run_[0-9]*"))):
    run = os.path.basename(d)
    s = pd.read_csv(os.path.join(d, "scores.tsv"), sep="\t", dtype=str)
    s.columns = [c.strip().lower() for c in s.columns]
    s = s.assign(pathway=s.pathway.str.strip()).set_index("pathway").score.astype(int)
    M[run] = s.reindex(M.index)
    r = int(run.split("_")[1])
    log.append(dict(run=run, shuffle_seed=man.shuffle_seed.get(r, ""),
                    **{f"n_score_{k}": int((s == k).sum()) for k in range(1, 6)}))
if M.drop(columns="term_name").isna().any().any():
    raise SystemExit("some GO IDs have no score; run validate_runs.py first")
os.makedirs(A.out_dir, exist_ok=True)
out = os.path.join(A.out_dir, f"score_matrix.{A.source}.tsv")
M.reset_index().to_csv(out, sep="\t", index=False)
pd.DataFrame(log).to_csv(os.path.join(src_dir, "run_log.tsv"), sep="\t", index=False)
print(f"wrote {out} ({len(M)} terms x {len(log)} runs)")
