#!/usr/bin/env python3
"""
Propagate chicken GO annotations to ancestor terms (true-path rule) and write a
gene-symbol -> GO term table in the same three-column format as the direct-annotation
file used by D1_map_GO_finch_chicken.py (gene_name, term_id, ontology; no header).

Propagation follows is_a and part_of edges (the standard choice; regulates edges
are not followed). Rows with a NOT qualifier are dropped before propagation, since a
NOT annotation must not be inherited by parents. The root term of each aspect is
excluded (it would contain every annotated gene).

Usage:
  python3 propagate_gaf.py --gaf gene_sets/CHICK-uniprot.gaf --obo gene_sets/go-basic.obo \
      --aspect P --out gene_sets/chicken_symbol_GO_BP.propagated.tsv

Requires: Python 3 only (go-basic.obo is parsed here; no goatools dependency).
"""
import argparse, gzip, sys
from collections import defaultdict


class Term:
    __slots__ = ("id", "name", "namespace", "is_a", "part_of", "is_obsolete", "replaced_by")
    def __init__(self):
        self.id = self.name = self.namespace = self.replaced_by = None
        self.is_a, self.part_of, self.is_obsolete = [], [], False


def read_obo(path):
    """Minimal parser for go-basic.obo: id, name, namespace, is_a, part_of, obsolete, replaced_by."""
    dag, t, in_term = {}, None, False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("["):
                if in_term and t.id: dag[t.id] = t
                in_term = line == "[Term]"; t = Term() if in_term else None
                continue
            if not in_term or not line or ":" not in line:
                continue
            key, val = line.split(":", 1); val = val.strip()
            if key == "id": t.id = val
            elif key == "name": t.name = val
            elif key == "namespace": t.namespace = val
            elif key == "is_a": t.is_a.append(val.split("!")[0].strip())
            elif key == "relationship":
                rel, tgt = val.split()[:2]
                if rel == "part_of": t.part_of.append(tgt)
            elif key == "is_obsolete": t.is_obsolete = val == "true"
            elif key == "replaced_by": t.replaced_by = val
        if in_term and t and t.id: dag[t.id] = t
    return dag

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--gaf", required=True, help="CHICK-uniprot.gaf (plain or .gz)")
ap.add_argument("--obo", required=True, help="go-basic.obo from the same GO release as the GAF")
ap.add_argument("--aspect", default="P", choices=["P", "F", "C"], help="GO aspect to keep (P = biological process)")
ap.add_argument("--out", required=True)
ap.add_argument("--direct-out", default=None,
                help="optional: also write the un-propagated table (same filters) for a like-for-like comparison")
A = ap.parse_args()

ASPECT_LABEL = {"P": "BP", "F": "MF", "C": "CC"}
ROOT = {"P": "GO:0008150", "F": "GO:0003674", "C": "GO:0005575"}

dag = read_obo(A.obo)
print(f"ontology terms read: {len(dag)}")

def ancestors(go_id, cache={}):
    """All ancestors of go_id via is_a and part_of, excluding go_id itself."""
    if go_id in cache:
        return cache[go_id]
    seen, stack = set(), [go_id]
    while stack:
        t = dag.get(stack.pop())
        if t is None:
            continue
        ups = set(t.is_a) | set(t.part_of)
        for u in ups:
            if u not in seen:
                seen.add(u); stack.append(u)
    cache[go_id] = seen
    return seen

opener = gzip.open if A.gaf.endswith(".gz") else open
direct = defaultdict(set)          # symbol -> {direct GO ids}
n_rows = n_not = n_obsolete = n_other_aspect = 0
with opener(A.gaf, "rt") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        n_rows += 1
        symbol, qualifier, go_id, aspect = f[2], f[3], f[4], f[8]
        if aspect != A.aspect:
            n_other_aspect += 1; continue
        if "NOT" in qualifier.split("|"):
            n_not += 1; continue
        term = dag.get(go_id)
        if term is None:
            n_obsolete += 1; continue
        if term.is_obsolete:
            # follow replaced_by if present, else drop
            rep = term.replaced_by
            if rep and rep in dag and not dag[rep].is_obsolete:
                go_id = rep
            else:
                n_obsolete += 1; continue
        direct[symbol].add(go_id)

label = ASPECT_LABEL[A.aspect]
root = ROOT[A.aspect]

propagated = {}
for sym, terms in direct.items():
    full = set(terms)
    for t in terms:
        full |= ancestors(t)
    full.discard(root)
    propagated[sym] = full

def write(table, path):
    n = 0
    with open(path, "w") as out:
        for sym in sorted(table):
            for go_id in sorted(table[sym]):
                out.write(f"{sym}\t{go_id}\t{label}\n"); n += 1
    return n

n_prop = write(propagated, A.out)
print(f"GAF rows: {n_rows}; other aspect: {n_other_aspect}; NOT dropped: {n_not}; obsolete/unknown dropped: {n_obsolete}")
print(f"symbols with {label} annotation: {len(direct)}")
print(f"direct gene-term pairs: {sum(len(v) for v in direct.values())}; "
      f"propagated gene-term pairs: {n_prop}")
print(f"wrote {A.out}")
if A.direct_out:
    for sym in direct: direct[sym].discard(root)
    n_dir = write(direct, A.direct_out)
    print(f"wrote {A.direct_out} ({n_dir} pairs)")
