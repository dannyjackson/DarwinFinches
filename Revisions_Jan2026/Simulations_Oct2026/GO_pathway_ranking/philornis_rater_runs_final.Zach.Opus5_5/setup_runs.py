#!/usr/bin/env python3
"""Set up 30 independent rater runs: shuffled inputs, filled prompts, seed manifest."""
import os, random

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")
N_RUNS = 30

PROMPT = """You are acting as an independent expert rater for a study of Darwin's finches and the
invasive avian vampire fly, Philornis downsi.

TASK
The file {INPUT} lists 1,008 Gene Ontology (GO) biological-process terms (one per line:
GO ID, a space, then the term name; the first line is a header). Score every term from
1 to 5 for how physiologically related it is to a potential host response to P. downsi.
The goal is an unbiased a priori list: if a term were later detected as under selection,
would we reasonably infer it is potentially associated with a response to the fly?

STEPS
1. Review the published literature on P. downsi and its effects on Galápagos land birds
   (Darwin's finches, mockingbirds and other passerines), using web search. Cover what the
   larvae do to hosts, and which host defenses (resistance, tolerance, behavioral) have been
   documented or proposed.
2. Base your scoring on that review.
3. Score every term.

SCALE
5 = directly captures a host mechanism documented in the Philornis literature
4 = core component of such a mechanism, or a strong a priori prediction from it
3 = plausible downstream or supporting process
2 = weak, indirect or speculative link
1 = no specific link

RULES
- Use only your literature review and general biological knowledge. Do not use memory,
  past conversations, or any information about results from genomic studies of these birds.
  Do not use transcriptomic, gene-expression or genomic studies of Galápagos birds (e.g., 
  Knutie et al. bioRxiv 2023.03.06.531275); if you encounter one, do not use its results.
- Judge each term on its own merits; do not aim for any particular distribution of scores.
- Score all 1,008 terms. Do not skip, merge or rename any.

INDEPENDENCE (required)
- Read only the file {INPUT}. Do not list, open or search any other file or directory on
  this machine, and do not use any project, document, memory or artifact tools. Other raters'
  work exists on this machine and you must not see it.
- Write only the two output files named below. If you write a helper script, put it in
  {RUNDIR} (never in a shared or temporary directory) and give it a name unique to this run.

OUTPUT
1. Write a brief summary of the mechanisms you identified and the key sources (with URLs)
   to {REVIEW}.
2. Write the scores to {SCORES} as a tab-separated file with the header "pathway<TAB>score"
   and exactly one row per GO ID in the input, using each GO ID exactly as given. Writing a
   short script that maps GO IDs to scores is fine, but every score must be your own judgment.
3. Reply with only one line: "done: <number of rows written>".
"""


def main():
    with open(SRC, encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f if l.strip()]
    header, terms = lines[0], lines[1:]
    assert len(terms) == 1008, len(terms)
    os.makedirs(RUNS, exist_ok=True)
    manifest = ["run\tseed"]
    for r in range(1, N_RUNS + 1):
        seed = 1000 + r
        rundir = os.path.join(RUNS, f"run_{r:03d}")
        os.makedirs(rundir, exist_ok=True)
        shuffled = terms[:]
        random.Random(seed).shuffle(shuffled)
        inp = os.path.join(rundir, "input_terms.txt")
        with open(inp, "w", encoding="utf-8") as f:
            f.write(header + "\n" + "\n".join(shuffled) + "\n")
        filled = (PROMPT.replace("{INPUT}", inp)
                        .replace("{REVIEW}", os.path.join(rundir, "review.md"))
                        .replace("{SCORES}", os.path.join(rundir, "scores.tsv"))
                        .replace("{RUNDIR}", rundir))
        with open(os.path.join(rundir, "prompt.txt"), "w", encoding="utf-8") as f:
            f.write(filled)
        manifest.append(f"run_{r:03d}\t{seed}")
    with open(os.path.join(RUNS, "run_manifest.tsv"), "w") as f:
        f.write("\n".join(manifest) + "\n")
    print("set up", N_RUNS, "runs")


if __name__ == "__main__":
    main()
