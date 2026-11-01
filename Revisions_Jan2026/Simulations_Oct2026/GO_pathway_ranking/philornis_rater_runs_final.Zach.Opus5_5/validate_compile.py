#!/usr/bin/env python3
"""Validate run outputs, and compile score_matrix.tsv / run_log.tsv.

Usage:
  python3 validate_compile.py validate run_001 run_002 ...
  python3 validate_compile.py compile
Token use and re-run flags are read from runs/run_status.tsv
(columns: run, tokens, rerun), which the coordinator maintains.
"""
import os, sys, csv

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")


def load_terms():
    with open(SRC, encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f if l.strip()][1:]
    return [tuple(l.split(" ", 1)) for l in lines]


def validate(run, ids):
    path = os.path.join(RUNS, run, "scores.tsv")
    if not os.path.exists(path):
        return False, "scores.tsv missing", None
    with open(path, encoding="utf-8") as f:
        rows = [l.rstrip("\r\n") for l in f if l.strip()]
    if not rows or rows[0].split("\t") != ["pathway", "score"]:
        return False, f"bad header: {rows[0] if rows else ''!r}", None
    data = rows[1:]
    if len(data) != 1008:
        return False, f"{len(data)} rows (need 1008)", None
    scores, problems = {}, []
    for l in data:
        parts = l.split("\t")
        if len(parts) != 2:
            problems.append(f"malformed row {l!r}"); continue
        gid, s = parts
        if gid in scores:
            problems.append(f"duplicate {gid}")
        if s not in {"1", "2", "3", "4", "5"}:
            problems.append(f"bad score {gid}={s!r}")
        scores[gid] = s
    extra = set(scores) - set(ids)
    missing = set(ids) - set(scores)
    if extra: problems.append(f"{len(extra)} extra IDs e.g. {sorted(extra)[:3]}")
    if missing: problems.append(f"{len(missing)} missing IDs e.g. {sorted(missing)[:3]}")
    if not os.path.exists(os.path.join(RUNS, run, "review.md")):
        problems.append("review.md missing")
    if problems:
        return False, "; ".join(problems[:5]), None
    return True, "ok", scores


def load_status():
    p = os.path.join(RUNS, "run_status.tsv")
    st = {}
    if os.path.exists(p):
        with open(p) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                st[row["run"]] = row
    return st


def compile_all():
    terms = load_terms()
    ids = [t[0] for t in terms]
    runs = [f"run_{r:03d}" for r in range(1, 31)]
    allscores = {}
    for run in runs:
        ok, msg, sc = validate(run, ids)
        if not ok:
            sys.exit(f"{run} fails validation: {msg}")
        allscores[run] = sc
    with open(os.path.join(BASE, "score_matrix.tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(["pathway", "term_name"] + runs) + "\n")
        for gid, name in terms:
            f.write("\t".join([gid, name] + [allscores[r][gid] for r in runs]) + "\n")
    seeds = {}
    with open(os.path.join(RUNS, "run_manifest.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            seeds[row["run"]] = row["seed"]
    st = load_status()
    with open(os.path.join(BASE, "run_log.tsv"), "w") as f:
        f.write("run\tseed\ttokens\trerun\tn_score1\tn_score2\tn_score3\tn_score4\tn_score5\n")
        for run in runs:
            vals = list(allscores[run].values())
            counts = [str(vals.count(str(k))) for k in range(1, 6)]
            s = st.get(run, {})
            f.write("\t".join([run, seeds[run], s.get("tokens", "NA"), s.get("rerun", "no")] + counts) + "\n")
    print("compiled score_matrix.tsv and run_log.tsv")


if __name__ == "__main__":
    if sys.argv[1] == "validate":
        ids = [t[0] for t in load_terms()]
        bad = 0
        for run in sys.argv[2:]:
            ok, msg, _ = validate(run, ids)
            print(f"{run}\t{'PASS' if ok else 'FAIL'}\t{msg}")
            bad += not ok
        sys.exit(1 if bad else 0)
    elif sys.argv[1] == "compile":
        compile_all()
