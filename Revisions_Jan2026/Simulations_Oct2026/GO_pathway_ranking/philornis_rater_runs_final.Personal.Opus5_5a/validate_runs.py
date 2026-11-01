#!/usr/bin/env python3
"""Validate rater outputs and make a backup zip of runs/ (input_terms.txt excluded).

Usage:
  python3 validate_runs.py validate 1 2 3 4 5
  python3 validate_runs.py backup <label>
"""
import os
import sys
import zipfile

BASE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(BASE, "all_pathways.txt")
RUNS = os.path.join(BASE, "runs")
BACKUPS = os.path.join(BASE, "backups")


def expected_ids():
    with open(SRC) as fh:
        lines = fh.read().splitlines()[1:]
    return [l.split(" ", 1)[0] for l in lines if l.strip()]


def validate(n):
    path = os.path.join(RUNS, f"run_{n:03d}", "scores.tsv")
    errs = []
    if not os.path.exists(path):
        return ["scores.tsv missing"]
    with open(path) as fh:
        lines = [l.rstrip("\r\n") for l in fh]
    while lines and not lines[-1].strip():
        lines.pop()
    if not lines or lines[0] != "pathway\tscore":
        errs.append(f"bad header: {lines[0]!r}" if lines else "empty file")
    rows = lines[1:]
    if len(rows) != 1008:
        errs.append(f"{len(rows)} data rows (expected 1008)")
    exp = set(expected_ids())
    seen = {}
    bad_scores = 0
    for r in rows:
        parts = r.split("\t")
        if len(parts) != 2:
            errs.append(f"malformed row: {r!r}")
            continue
        gid, sc = parts
        seen[gid] = seen.get(gid, 0) + 1
        if sc not in {"1", "2", "3", "4", "5"}:
            bad_scores += 1
    dups = [g for g, c in seen.items() if c > 1]
    missing = exp - set(seen)
    extra = set(seen) - exp
    if dups:
        errs.append(f"{len(dups)} duplicated IDs")
    if missing:
        errs.append(f"{len(missing)} missing IDs")
    if extra:
        errs.append(f"{len(extra)} extra IDs: {sorted(extra)[:5]}")
    if bad_scores:
        errs.append(f"{bad_scores} non-integer or out-of-range scores")
    if not os.path.exists(os.path.join(RUNS, f"run_{n:03d}", "review.md")):
        errs.append("review.md missing")
    return errs


def backup(label):
    os.makedirs(BACKUPS, exist_ok=True)
    out = os.path.join(BACKUPS, f"runs_backup_{label}.zip")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(RUNS):
            for f in sorted(files):
                if f == "input_terms.txt":
                    continue
                full = os.path.join(root, f)
                z.write(full, os.path.relpath(full, BASE))
    print(out)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "validate":
        for a in sys.argv[2:]:
            n = int(a)
            e = validate(n)
            print(f"run_{n:03d}\t{'PASS' if not e else 'FAIL: ' + '; '.join(e)}")
    elif cmd == "backup":
        backup(sys.argv[2])
