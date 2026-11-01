#!/usr/bin/env python3
"""Setup, validation and compilation for the 30 independent GO-term scoring runs.

Usage:
  setup_and_compile.py setup    <all_pathways.txt> <runs_dir> <prompt_template> <user_facing_runs_dir>
  setup_and_compile.py validate <all_pathways.txt> <runs_dir> [run_numbers...]
  setup_and_compile.py compile  <all_pathways.txt> <runs_dir> <out_dir>

<user_facing_runs_dir> is the path to runs/ as the sub-agents' file tools see it
(it may differ from the path used by this script if the script runs inside a VM mount).
"""
import csv
import os
import random
import sys
from collections import Counter

N_RUNS = 30
SEED_BASE = 1000


def read_pathways(path):
    with open(path, encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    header, body = lines[0], lines[1:]
    ids = [ln.split(" ", 1)[0] for ln in body]
    names = [ln.split(" ", 1)[1] if " " in ln else "" for ln in body]
    assert len(set(ids)) == len(ids), "duplicate GO IDs in input"
    return header, body, ids, names


def setup(all_pathways, runs_dir, template_path, user_runs_dir):
    header, body, ids, _ = read_pathways(all_pathways)
    with open(template_path, encoding="utf-8") as fh:
        template = fh.read()
    os.makedirs(runs_dir, exist_ok=True)
    manifest = []
    for n in range(1, N_RUNS + 1):
        name = f"run_{n:03d}"
        seed = SEED_BASE + n
        rundir = os.path.join(runs_dir, name)
        os.makedirs(rundir, exist_ok=True)
        terms = list(body)
        random.Random(seed).shuffle(terms)
        with open(os.path.join(rundir, "input_terms.txt"), "w", encoding="utf-8") as fh:
            fh.write(header + "\n")
            fh.write("\n".join(terms) + "\n")
        urd = os.path.join(user_runs_dir, name)
        prompt = (template.replace("{INPUT}", os.path.join(urd, "input_terms.txt"))
                          .replace("{REVIEW}", os.path.join(urd, "review.md"))
                          .replace("{SCORES}", os.path.join(urd, "scores.tsv"))
                          .replace("{RUNDIR}", urd))
        with open(os.path.join(rundir, "prompt.txt"), "w", encoding="utf-8") as fh:
            fh.write(prompt)
        manifest.append((name, n, seed, len(terms)))
    with open(os.path.join(runs_dir, "run_manifest.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "run_number", "shuffle_seed", "n_terms"])
        w.writerows(manifest)
    print(f"set up {N_RUNS} runs in {runs_dir}")


def validate_one(ids, scores_path):
    """Return (ok, message, scores_dict)."""
    expected = set(ids)
    if not os.path.exists(scores_path):
        return False, "scores.tsv missing", None
    with open(scores_path, encoding="utf-8") as fh:
        rows = [ln.rstrip("\r\n") for ln in fh]
    # tolerate trailing blank line(s) only
    while rows and rows[-1] == "":
        rows.pop()
    if not rows or rows[0].split("\t") != ["pathway", "score"]:
        return False, f"bad header: {rows[0] if rows else '<empty>'}", None
    data = rows[1:]
    problems = []
    if len(data) != len(ids):
        problems.append(f"{len(data)} rows (expected {len(ids)})")
    seen = Counter()
    scores = {}
    for i, ln in enumerate(data, start=2):
        parts = ln.split("\t")
        if len(parts) != 2:
            problems.append(f"line {i}: {len(parts)} fields")
            continue
        gid, sc = parts
        seen[gid] += 1
        if not (sc.isdigit() and 1 <= int(sc) <= 5):
            problems.append(f"line {i}: non-integer or out-of-range score '{sc}'")
            continue
        scores[gid] = int(sc)
    dup = [g for g, c in seen.items() if c > 1]
    extra = sorted(set(seen) - expected)
    missing = sorted(expected - set(seen))
    if dup:
        problems.append(f"{len(dup)} duplicated IDs (e.g. {dup[:3]})")
    if extra:
        problems.append(f"{len(extra)} extra IDs (e.g. {extra[:3]})")
    if missing:
        problems.append(f"{len(missing)} missing IDs (e.g. {missing[:3]})")
    if problems:
        return False, "; ".join(problems), None
    return True, f"OK ({len(data)} rows)", scores


def validate(all_pathways, runs_dir, run_numbers):
    _, _, ids, _ = read_pathways(all_pathways)
    nums = [int(x) for x in run_numbers] if run_numbers else range(1, N_RUNS + 1)
    all_ok = True
    for n in nums:
        name = f"run_{n:03d}"
        ok, msg, _ = validate_one(ids, os.path.join(runs_dir, name, "scores.tsv"))
        all_ok &= ok
        print(f"{name}\t{'PASS' if ok else 'FAIL'}\t{msg}")
    sys.exit(0 if all_ok else 1)


def compile_results(all_pathways, runs_dir, out_dir):
    _, _, ids, names = read_pathways(all_pathways)
    manifest = {}
    with open(os.path.join(runs_dir, "run_manifest.tsv"), encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            manifest[row["run"]] = row
    # tokens / rerun info recorded by the coordinator in runs/run_status.tsv
    status = {}
    status_path = os.path.join(runs_dir, "run_status.tsv")
    if os.path.exists(status_path):
        with open(status_path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                status[row["run"]] = row
    run_names = [f"run_{n:03d}" for n in range(1, N_RUNS + 1)]
    per_run = {}
    for name in run_names:
        ok, msg, scores = validate_one(ids, os.path.join(runs_dir, name, "scores.tsv"))
        if not ok:
            sys.exit(f"{name} fails validation: {msg}")
        per_run[name] = scores
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "score_matrix.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pathway", "term_name"] + run_names)
        for gid, nm in zip(ids, names):
            w.writerow([gid, nm] + [per_run[r][gid] for r in run_names])
    with open(os.path.join(out_dir, "run_log.tsv"), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["run", "run_number", "shuffle_seed", "token_use", "rerun",
                    "n_score_1", "n_score_2", "n_score_3", "n_score_4", "n_score_5"])
        for name in run_names:
            c = Counter(per_run[name].values())
            st = status.get(name, {})
            w.writerow([name, manifest[name]["run_number"], manifest[name]["shuffle_seed"],
                        st.get("token_use", "NA"), st.get("rerun", "NA")]
                       + [c.get(s, 0) for s in range(1, 6)])
    print(f"wrote score_matrix.tsv ({len(ids)} x {N_RUNS}) and run_log.tsv to {out_dir}")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "setup":
        setup(*sys.argv[2:6])
    elif cmd == "validate":
        validate(sys.argv[2], sys.argv[3], sys.argv[4:])
    elif cmd == "compile":
        compile_results(*sys.argv[2:5])
    else:
        sys.exit(__doc__)
