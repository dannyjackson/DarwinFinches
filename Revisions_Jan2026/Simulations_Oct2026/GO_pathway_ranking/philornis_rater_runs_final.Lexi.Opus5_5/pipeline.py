#!/usr/bin/env python3
"""Setup / validate / backup / compile script for the 30-run P. downsi GO-term rating study.

Usage:
  python3 pipeline.py setup
  python3 pipeline.py validate RUN [RUN ...]
  python3 pipeline.py backup NAME.zip
  python3 pipeline.py compile
Run-level metadata (token use, re-runs) is kept in runs/run_status.tsv.
"""
import csv, os, random, sys, zipfile, re
from collections import Counter

BASE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(BASE, "runs")
PATHWAYS = "/mnt/user-data/uploads/607fdcaf-all_pathways.txt"
INSTRUCTIONS = "/mnt/user-data/uploads/3eb75a1e-attachment.txt"
N_RUNS = 30
STATUS = os.path.join(RUNS, "run_status.tsv")


def run_dir(n):
    return os.path.join(RUNS, f"run_{n:03d}")


def read_pathways():
    with open(PATHWAYS, encoding="utf-8") as f:
        lines = [l.rstrip("\n").rstrip("\r") for l in f]
    header, body = lines[0], [l for l in lines[1:] if l.strip()]
    terms = []
    for l in body:
        gid, name = l.split(" ", 1)
        terms.append((gid, name))
    return header, body, terms


def rater_template():
    """Extract the rater prompt verbatim from the instructions file (between the --- lines)."""
    with open(INSTRUCTIONS, encoding="utf-8") as f:
        text = f.read()
    start = text.index("RATER PROMPT")
    seg = text[start:]
    a = seg.index("---\n") + 4
    b = seg.index("\n---", a)
    return seg[a:b] + "\n"


def setup():
    header, body, terms = read_pathways()
    assert len(body) == 1008, len(body)
    tmpl = rater_template()
    for ph in ("{INPUT}", "{REVIEW}", "{SCORES}", "{RUNDIR}"):
        assert ph in tmpl, ph
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "run_manifest.tsv"), "w") as m:
        m.write("run\tseed\trun_dir\n")
        for n in range(1, N_RUNS + 1):
            d = run_dir(n)
            os.makedirs(d, exist_ok=True)
            seed = 1000 + n
            lines = list(body)
            random.Random(seed).shuffle(lines)
            inp = os.path.join(d, "input_terms.txt")
            with open(inp, "w", encoding="utf-8") as f:
                f.write(header + "\n" + "\n".join(lines) + "\n")
            filled = (tmpl.replace("{INPUT}", inp)
                          .replace("{REVIEW}", os.path.join(d, "review.md"))
                          .replace("{SCORES}", os.path.join(d, "scores.tsv"))
                          .replace("{RUNDIR}", d))
            with open(os.path.join(d, "prompt.txt"), "w", encoding="utf-8") as f:
                f.write(filled)
            m.write(f"run_{n:03d}\t{seed}\t{d}\n")
    if not os.path.exists(STATUS):
        with open(STATUS, "w") as f:
            f.write("run\tattempts\ttotal_tokens\tnote\n")
    print("setup ok")


def validate(n):
    _, _, terms = read_pathways()
    ids = [g for g, _ in terms]
    idset = set(ids)
    p = os.path.join(run_dir(n), "scores.tsv")
    errs = []
    if not os.path.exists(p):
        return [f"missing {p}"]
    with open(p, encoding="utf-8") as f:
        rows = [l.rstrip("\n").rstrip("\r") for l in f if l.strip()]
    if not rows or rows[0].split("\t") != ["pathway", "score"]:
        errs.append(f"bad header: {rows[0] if rows else None!r}")
    data = rows[1:]
    if len(data) != 1008:
        errs.append(f"{len(data)} rows (expected 1008)")
    seen = Counter()
    for r in data:
        parts = r.split("\t")
        if len(parts) != 2:
            errs.append(f"bad row {r!r}"); continue
        g, s = parts
        seen[g] += 1
        if g not in idset:
            errs.append(f"extra id {g!r}")
        if not re.fullmatch(r"[1-5]", s.strip()):
            errs.append(f"bad score {g}:{s!r}")
    dups = [g for g, c in seen.items() if c > 1]
    missing = [g for g in ids if g not in seen]
    if dups: errs.append(f"duplicates: {dups[:10]}")
    if missing: errs.append(f"missing {len(missing)} ids: {missing[:10]}")
    return errs


def load_scores(n):
    with open(os.path.join(run_dir(n), "scores.tsv"), encoding="utf-8") as f:
        rows = [l.rstrip("\n").rstrip("\r").split("\t") for l in f if l.strip()][1:]
    return {g: int(s.strip()) for g, s in rows}


def backup(name):
    out = os.path.join(BASE, name)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(RUNS):
            for fn in sorted(files):
                if fn == "input_terms.txt":
                    continue
                full = os.path.join(root, fn)
                z.write(full, os.path.relpath(full, BASE))
    print(out)


def read_status():
    st = {}
    with open(STATUS) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st[r["run"]] = r
    return st


def compile_all():
    _, _, terms = read_pathways()
    for n in range(1, N_RUNS + 1):
        e = validate(n)
        assert not e, (n, e)
    scores = {n: load_scores(n) for n in range(1, N_RUNS + 1)}
    with open(os.path.join(BASE, "score_matrix.tsv"), "w", encoding="utf-8") as f:
        f.write("\t".join(["pathway", "term_name"] + [f"run_{n:03d}" for n in range(1, N_RUNS + 1)]) + "\n")
        for g, name in terms:
            f.write("\t".join([g, name] + [str(scores[n][g]) for n in range(1, N_RUNS + 1)]) + "\n")
    st = read_status()
    with open(os.path.join(BASE, "run_log.tsv"), "w") as f:
        f.write("run\tseed\ttotal_tokens\trerun\tattempts\tn_score1\tn_score2\tn_score3\tn_score4\tn_score5\tnote\n")
        for n in range(1, N_RUNS + 1):
            key = f"run_{n:03d}"
            s = st.get(key, {})
            c = Counter(scores[n].values())
            att = int(s.get("attempts", 1) or 1)
            f.write("\t".join([key, str(1000 + n), s.get("total_tokens", ""),
                               "yes" if att > 1 else "no", str(att)]
                              + [str(c.get(k, 0)) for k in range(1, 6)] + [s.get("note", "")]) + "\n")
    name = os.path.join(BASE, "final_results.zip")
    with zipfile.ZipFile(name, "w", zipfile.ZIP_DEFLATED) as z:
        for root, _, files in os.walk(RUNS):
            for fn in sorted(files):
                full = os.path.join(root, fn)
                z.write(full, os.path.relpath(full, BASE))
        for fn in ("score_matrix.tsv", "run_log.tsv", "pipeline.py"):
            z.write(os.path.join(BASE, fn), fn)
    print(name)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "setup":
        setup()
    elif cmd == "validate":
        for a in sys.argv[2:]:
            e = validate(int(a))
            print(f"run_{int(a):03d}: {'PASS' if not e else 'FAIL ' + '; '.join(e)}")
    elif cmd == "backup":
        backup(sys.argv[2])
    elif cmd == "compile":
        compile_all()
