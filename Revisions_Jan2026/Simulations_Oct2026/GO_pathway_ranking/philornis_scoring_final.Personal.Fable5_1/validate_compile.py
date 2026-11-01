#!/usr/bin/env python3
"""Validate scores.tsv files and (optionally) compile score_matrix.tsv / run_log.tsv.

Usage:
  validate_compile.py validate run_001 run_002 ...   # validate listed runs (or all if none given)
  validate_compile.py compile                        # build score_matrix.tsv and run_log.tsv
"""
import csv
import json
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")
STATUS = os.path.join(RUNS, "run_status.json")   # token use / re-run flags, maintained by coordinator
N_RUNS = 30

with open(SRC) as f:
    lines = [l.rstrip("\n") for l in f if l.strip()]
ORDER = [l.split(" ", 1)[0] for l in lines[1:]]
NAMES = {l.split(" ", 1)[0]: l.split(" ", 1)[1] for l in lines[1:]}
EXPECTED = set(ORDER)
assert len(ORDER) == 1008 and len(EXPECTED) == 1008


def validate(run):
    path = os.path.join(RUNS, run, "scores.tsv")
    problems = []
    if not os.path.exists(path):
        return ["scores.tsv missing"], None
    with open(path) as f:
        rows = [l.rstrip("\n") for l in f if l.strip()]
    if not rows:
        return ["empty file"], None
    hdr = rows[0].split("\t")
    if hdr != ["pathway", "score"]:
        problems.append(f"bad header: {rows[0]!r}")
    data = rows[1:]
    if len(data) != 1008:
        problems.append(f"{len(data)} data rows (expected 1008)")
    scores = {}
    seen = []
    for r in data:
        parts = r.split("\t")
        if len(parts) != 2:
            problems.append(f"bad row: {r!r}")
            continue
        gid, sc = parts
        seen.append(gid)
        try:
            v = int(sc)
            if not 1 <= v <= 5:
                problems.append(f"{gid}: score {sc} out of range")
        except ValueError:
            problems.append(f"{gid}: non-integer score {sc!r}")
            continue
        scores[gid] = v
    dup = {g for g in seen if seen.count(g) > 1}
    if dup:
        problems.append(f"duplicate IDs: {sorted(dup)[:10]}")
    missing = EXPECTED - set(seen)
    extra = set(seen) - EXPECTED
    if missing:
        problems.append(f"missing {len(missing)} IDs, e.g. {sorted(missing)[:5]}")
    if extra:
        problems.append(f"{len(extra)} extra IDs, e.g. {sorted(extra)[:5]}")
    if not os.path.exists(os.path.join(RUNS, run, "review.md")):
        problems.append("review.md missing")
    return problems, scores


def load_status():
    if os.path.exists(STATUS):
        with open(STATUS) as f:
            return json.load(f)
    return {}


def compile_all():
    status = load_status()
    manifest = {}
    with open(os.path.join(RUNS, "run_manifest.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            manifest[row["run"]] = row["seed"]
    runs = [f"run_{i:03d}" for i in range(1, N_RUNS + 1)]
    all_scores = {}
    for run in runs:
        problems, scores = validate(run)
        if problems:
            sys.exit(f"{run} fails validation: {problems}")
        all_scores[run] = scores
    with open(os.path.join(BASE, "score_matrix.tsv"), "w") as f:
        f.write("pathway\tterm_name\t" + "\t".join(runs) + "\n")
        for gid in ORDER:
            f.write(gid + "\t" + NAMES[gid] + "\t" + "\t".join(str(all_scores[r][gid]) for r in runs) + "\n")
    with open(os.path.join(BASE, "run_log.tsv"), "w") as f:
        f.write("run\tseed\ttokens_total\ttokens_note\trerun\tn_1\tn_2\tn_3\tn_4\tn_5\n")
        for run in runs:
            st = status.get(run, {})
            counts = [sum(1 for v in all_scores[run].values() if v == k) for k in range(1, 6)]
            f.write("\t".join([run, manifest[run], str(st.get("tokens_total", "")),
                               str(st.get("tokens_note", "")),
                               "yes" if st.get("rerun") else "no"] + [str(c) for c in counts]) + "\n")
    print("wrote score_matrix.tsv and run_log.tsv")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "compile":
        compile_all()
    else:
        targets = sys.argv[2:] if len(sys.argv) > 1 and sys.argv[1] == "validate" else sys.argv[1:]
        if not targets:
            targets = [f"run_{i:03d}" for i in range(1, N_RUNS + 1)]
        ok = True
        for run in targets:
            problems, scores = validate(run)
            if problems:
                ok = False
                print(f"FAIL {run}: " + "; ".join(problems))
            else:
                counts = [sum(1 for v in scores.values() if v == k) for k in range(1, 6)]
                print(f"PASS {run}: 1008 rows, counts 1-5 = {counts}")
        sys.exit(0 if ok else 1)
