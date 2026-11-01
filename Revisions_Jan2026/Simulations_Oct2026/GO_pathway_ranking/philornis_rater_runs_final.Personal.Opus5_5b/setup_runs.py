#!/usr/bin/env python3
"""Set up 30 independent rater runs and compile their scores.

Usage:
  python3 setup_runs.py setup              # create runs/run_001..run_030, inputs, prompts, manifest
  python3 setup_runs.py validate N [N...]  # validate scores.tsv for given run numbers
  python3 setup_runs.py compile            # build score_matrix.tsv and run_log.tsv
"""
import csv
import os
import random
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")
N_RUNS = 30
N_TERMS = 1008

RATER_PROMPT = """You are acting as an independent expert rater for a study of Darwin's finches and the
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


def read_terms():
    with open(SRC, encoding="utf-8") as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    header, body = lines[0], lines[1:]
    assert len(body) == N_TERMS, f"expected {N_TERMS} terms, got {len(body)}"
    ids = [l.split(" ", 1)[0] for l in body]
    assert len(set(ids)) == N_TERMS, "duplicate GO IDs in source"
    return header, body


def run_dir(n):
    return os.path.join(RUNS, f"run_{n:03d}")


def setup():
    header, body = read_terms()
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "run_manifest.tsv"), "w", encoding="utf-8") as man:
        man.write("run\tseed\n")
        for n in range(1, N_RUNS + 1):
            d = run_dir(n)
            os.makedirs(d, exist_ok=True)
            seed = 1000 + n
            terms = list(body)
            random.Random(seed).shuffle(terms)
            with open(os.path.join(d, "input_terms.txt"), "w", encoding="utf-8") as fh:
                fh.write(header + "\n" + "\n".join(terms) + "\n")
            prompt = (RATER_PROMPT
                      .replace("{INPUT}", os.path.join(d, "input_terms.txt"))
                      .replace("{REVIEW}", os.path.join(d, "review.md"))
                      .replace("{SCORES}", os.path.join(d, "scores.tsv"))
                      .replace("{RUNDIR}", d))
            with open(os.path.join(d, "prompt.txt"), "w", encoding="utf-8") as fh:
                fh.write(prompt)
            man.write(f"run_{n:03d}\t{seed}\n")
    print(f"set up {N_RUNS} runs in {RUNS}")


def validate_one(n, expected):
    path = os.path.join(run_dir(n), "scores.tsv")
    if not os.path.exists(path):
        return False, "scores.tsv missing"
    with open(path, encoding="utf-8") as fh:
        rows = [l.rstrip("\r\n") for l in fh if l.strip()]
    if not rows or rows[0].split("\t") != ["pathway", "score"]:
        return False, f"bad header: {rows[0] if rows else '(empty)'}"
    data = rows[1:]
    problems = []
    if len(data) != N_TERMS:
        problems.append(f"{len(data)} rows (expected {N_TERMS})")
    seen = {}
    for i, r in enumerate(data, 2):
        parts = r.split("\t")
        if len(parts) != 2:
            problems.append(f"line {i}: {len(parts)} fields")
            continue
        gid, s = parts
        seen[gid] = seen.get(gid, 0) + 1
        if s not in {"1", "2", "3", "4", "5"}:
            problems.append(f"line {i}: score {s!r}")
    dups = [g for g, c in seen.items() if c > 1]
    missing = expected - set(seen)
    extra = set(seen) - expected
    if dups:
        problems.append(f"{len(dups)} duplicated IDs")
    if missing:
        problems.append(f"{len(missing)} missing IDs")
    if extra:
        problems.append(f"{len(extra)} extra IDs")
    return (not problems), "; ".join(problems[:6]) or "ok"


def validate(nums):
    _, body = read_terms()
    expected = {l.split(" ", 1)[0] for l in body}
    all_ok = True
    for n in nums:
        ok, msg = validate_one(n, expected)
        all_ok &= ok
        print(f"run_{n:03d}\t{'PASS' if ok else 'FAIL'}\t{msg}")
    return all_ok


def compile_all():
    _, body = read_terms()
    order = [l.split(" ", 1) for l in body]
    scores = {}
    for n in range(1, N_RUNS + 1):
        with open(os.path.join(run_dir(n), "scores.tsv"), encoding="utf-8") as fh:
            rdr = csv.reader(fh, delimiter="\t")
            next(rdr)
            scores[n] = {r[0]: r[1] for r in rdr if r}
    with open(os.path.join(BASE, "score_matrix.tsv"), "w", encoding="utf-8") as out:
        out.write("\t".join(["pathway", "term_name"] + [f"run_{n:03d}" for n in range(1, N_RUNS + 1)]) + "\n")
        for gid, name in order:
            out.write("\t".join([gid, name] + [scores[n][gid] for n in range(1, N_RUNS + 1)]) + "\n")
    # token use and re-run flags come from runs/run_status.tsv (kept by the coordinator)
    status = {}
    sp = os.path.join(RUNS, "run_status.tsv")
    if os.path.exists(sp):
        with open(sp, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                status[r["run"]] = r
    with open(os.path.join(BASE, "run_log.tsv"), "w", encoding="utf-8") as out:
        out.write("run\tseed\ttokens\trerun\tn_score1\tn_score2\tn_score3\tn_score4\tn_score5\n")
        for n in range(1, N_RUNS + 1):
            key = f"run_{n:03d}"
            st = status.get(key, {})
            counts = [sum(1 for v in scores[n].values() if v == str(k)) for k in range(1, 6)]
            out.write("\t".join([key, str(1000 + n), st.get("tokens", "NA"), st.get("rerun", "no")]
                                + [str(c) for c in counts]) + "\n")
    print("wrote score_matrix.tsv and run_log.tsv")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "setup":
        setup()
    elif cmd == "validate":
        sys.exit(0 if validate([int(x) for x in sys.argv[2:]]) else 1)
    elif cmd == "compile":
        compile_all()
    else:
        print(__doc__)
        sys.exit(2)
