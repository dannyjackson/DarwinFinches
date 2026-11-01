#!/usr/bin/env python3
"""Set up 30 independent rater runs.

Creates runs/run_001 ... runs/run_030, each containing:
  - input_terms.txt : header + all 1,008 terms shuffled with random.Random(1000 + run_number)
  - prompt.txt      : the rater prompt with {INPUT}, {REVIEW}, {SCORES}, {RUNDIR} filled in
and runs/run_manifest.tsv logging the seed for each run.
"""
import os
import random
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")
TEMPLATE = os.path.join(BASE, "rater_prompt_template.txt")
N_RUNS = 30

with open(SRC) as f:
    lines = [l.rstrip("\n") for l in f if l.strip()]
header, terms = lines[0], lines[1:]
assert len(terms) == 1008, len(terms)

with open(TEMPLATE) as f:
    template = f.read()

os.makedirs(RUNS, exist_ok=True)
manifest = ["run\tseed\tinput_file"]
for i in range(1, N_RUNS + 1):
    run = f"run_{i:03d}"
    rundir = os.path.join(RUNS, run)
    os.makedirs(rundir, exist_ok=True)
    seed = 1000 + i
    shuffled = list(terms)
    random.Random(seed).shuffle(shuffled)
    inp = os.path.join(rundir, "input_terms.txt")
    with open(inp, "w") as f:
        f.write(header + "\n" + "\n".join(shuffled) + "\n")
    review = os.path.join(rundir, "review.md")
    scores = os.path.join(rundir, "scores.tsv")
    prompt = (template.replace("{INPUT}", inp)
                      .replace("{REVIEW}", review)
                      .replace("{SCORES}", scores)
                      .replace("{RUNDIR}", rundir))
    with open(os.path.join(rundir, "prompt.txt"), "w") as f:
        f.write(prompt)
    manifest.append(f"{run}\t{seed}\t{inp}")

with open(os.path.join(RUNS, "run_manifest.tsv"), "w") as f:
    f.write("\n".join(manifest) + "\n")
print(f"set up {N_RUNS} runs in {RUNS}")
