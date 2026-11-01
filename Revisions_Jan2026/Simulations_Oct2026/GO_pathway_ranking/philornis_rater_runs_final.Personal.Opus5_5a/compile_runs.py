#!/usr/bin/env python3
"""Compile the 30 validated rater runs.

Inputs:  all_pathways.txt, runs/run_manifest.tsv, runs/run_status.tsv
         (run, total_tokens, rerun), runs/run_NNN/scores.tsv
Outputs: score_matrix.tsv  (pathway, term_name, run_001..run_030; all_pathways.txt order)
         run_log.tsv       (run, seed, total_tokens, rerun, n_score_1..n_score_5)
No analysis or interpretation of scores is done here.
"""
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(BASE, "runs")
N_RUNS = 30


def read_tsv(path):
    with open(path) as fh:
        rows = [l.rstrip("\r\n").split("\t") for l in fh if l.strip()]
    return rows[0], rows[1:]


def main():
    with open(os.path.join(BASE, "all_pathways.txt")) as fh:
        lines = [l for l in fh.read().splitlines()[1:] if l.strip()]
    terms = [l.split(" ", 1) for l in lines]
    ids = [t[0] for t in terms]

    _, man = read_tsv(os.path.join(RUNS, "run_manifest.tsv"))
    seeds = {r[0]: r[1] for r in man}
    _, st = read_tsv(os.path.join(RUNS, "run_status.tsv"))
    status = {r[0]: r[1:] for r in st}

    runs = [f"run_{n:03d}" for n in range(1, N_RUNS + 1)]
    scores = {}
    for run in runs:
        _, rows = read_tsv(os.path.join(RUNS, run, "scores.tsv"))
        d = {r[0]: r[1] for r in rows}
        if set(d) != set(ids) or len(rows) != len(ids):
            sys.exit(f"{run}: IDs do not match all_pathways.txt")
        scores[run] = d

    with open(os.path.join(BASE, "score_matrix.tsv"), "w") as fh:
        fh.write("\t".join(["pathway", "term_name"] + runs) + "\n")
        for gid, name in terms:
            fh.write("\t".join([gid, name] + [scores[r][gid] for r in runs]) + "\n")

    with open(os.path.join(BASE, "run_log.tsv"), "w") as fh:
        fh.write("\t".join(["run", "seed", "total_tokens", "rerun"]
                           + [f"n_score_{k}" for k in range(1, 6)]) + "\n")
        for run in runs:
            vals = list(scores[run].values())
            counts = [str(vals.count(str(k))) for k in range(1, 6)]
            tok, rerun = status.get(run, ["NA", "NA"])
            fh.write("\t".join([run[4:], seeds[run], tok, rerun] + counts) + "\n")
    print("wrote score_matrix.tsv and run_log.tsv")


if __name__ == "__main__":
    main()
