#!/usr/bin/env python3
"""
Set up the run folders for one source (one model on one account).

For each run r = 1..N this writes runs/<source>/run_<rrr>/:
  input_terms.txt  the header plus all terms in a random order (Python random.Random(1000 + r).shuffle())
  prompt.txt       the rater prompt with {INPUT}, {REVIEW}, {SCORES} and {RUNDIR} filled in for this run
and runs/<source>/run_manifest.tsv with the shuffle seed of each run.

This is STEP 1 of prompts/coordinating_prompt.txt. A coordinating conversation can run it, or it can be run by hand
before giving each prompt.txt to a fresh, independent rater (one new conversation or sub-agent per run).

Usage:
  export PROJDIR=/path/to/project
  python3 ${PROJDIR}/LLM_evaluation/scripts/setup_runs.py --source informed.Opus5_5.run_a --n-runs 30
"""
import argparse, os, random

ROOT = os.path.join(os.environ["PROJDIR"], "LLM_evaluation") if "PROJDIR" in os.environ else \
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--source", required=True, help="name for this model x account, e.g. naive.Fable5_1.Rebecca")
ap.add_argument("--n-runs", type=int, default=30)
ap.add_argument("--terms", default=os.path.join(ROOT, "runs", "all_pathways.txt"))
ap.add_argument("--prompt", default=os.path.join(ROOT, "prompts", "rater_prompt.txt"))
ap.add_argument("--runs-dir", default=os.path.join(ROOT, "runs"))
A = ap.parse_args()

lines = open(A.terms).read().splitlines()
header, terms = lines[0], [l for l in lines[1:] if l.strip()]
template = open(A.prompt).read()
src_dir = os.path.join(A.runs_dir, A.source)
os.makedirs(src_dir, exist_ok=True)
with open(os.path.join(src_dir, "run_manifest.tsv"), "w") as man:
    man.write("run\tshuffle_seed\n")
    for r in range(1, A.n_runs + 1):
        d = os.path.abspath(os.path.join(src_dir, f"run_{r:03d}"))
        os.makedirs(d, exist_ok=True)
        t = terms[:]
        random.Random(1000 + r).shuffle(t)
        inp = os.path.join(d, "input_terms.txt")
        open(inp, "w").write("\n".join([header] + t) + "\n")
        filled = (template.replace("{INPUT}", inp).replace("{REVIEW}", os.path.join(d, "review.md"))
                  .replace("{SCORES}", os.path.join(d, "scores.tsv")).replace("{RUNDIR}", d))
        open(os.path.join(d, "prompt.txt"), "w").write(filled)
        man.write(f"{r}\t{1000 + r}\n")
print(f"set up {A.n_runs} runs ({len(terms)} terms each) in {src_dir}")
