#!/usr/bin/env python3
"""
Compute temporal statistics from a SLiM run that started from an
msprime-built population (make_initial_state.py).

No recapitation is needed: the starting genealogy was complete. The script
checks this and recapitates only if some tree still has multiple roots.

Optionally overlays neutral mutations on the branches of the forward-simulated
period (the msprime overlay in make_initial_state.py covers everything older).

Output TSV keeps the same leading columns as recap_fst_tajd.all_parameter_space.py
so existing plotting scripts keep working; new columns are appended at the end.
"""
import argparse
import os

import msprime
import numpy as np
import pyslim
import tskit

# msprime >= 1.4 renames the SLiM model; older releases only have SLiMMutationModel.
# These mutations are only used for statistics here, never loaded back into SLiM.
SLIM_MUT_MODEL = getattr(msprime, "SLiMv5MutationModel", None) or msprime.SLiMMutationModel


def load_ids(path):
    with open(path) as f:
        return {int(x) for x in f.read().split()}


def sample_nodes_from_pedigrees(ts, pedigree_ids):
    nodes = []
    for ind in ts.individuals():
        if ind.metadata["pedigree_id"] in pedigree_ids:
            nodes.extend(n for n in ind.nodes if n != tskit.NULL)
    if not nodes:
        raise RuntimeError(f"No nodes found for {len(pedigree_ids)} pedigree IDs")
    return np.array(nodes, dtype=np.int32)


def _f(x):
    """tskit returns a scalar or a length-1 array depending on version."""
    return float(np.asarray(x).ravel()[0])


def hudson_fst(ts, a, b):
    pi_a = _f(ts.diversity(sample_sets=[a]))
    pi_b = _f(ts.diversity(sample_sets=[b]))
    dxy = _f(ts.divergence(sample_sets=[a, b]))
    fst = np.nan if dxy == 0 else 1.0 - 0.5 * (pi_a + pi_b) / dxy
    return fst, pi_a, pi_b, dxy


def tajimas_d(ts, nodes):
    sub = ts.simplify(samples=nodes, keep_unary=True)
    return _f(sub.Tajimas_D())


def site_freq(ts, pos, nodes):
    """Derived-allele frequency of the site at `pos` among `nodes` (0 if absent)."""
    sites = [s for s in ts.sites() if int(s.position) == pos]
    if not sites:
        return np.nan
    site = sites[0]
    var = next(ts.variants(samples=nodes, left=site.position, right=site.position + 1))
    return float(np.mean(var.genotypes != 0))


def max_roots(ts):
    return max(t.num_roots for t in ts.trees())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--trees", required=True)
    ap.add_argument("--t1_ids", required=True)
    ap.add_argument("--t2_ids", required=True)
    ap.add_argument("--site_file", required=True, help="from make_initial_state.py")
    ap.add_argument("--Ne", type=float, required=True)
    ap.add_argument("--mu", type=float, default=2.04e-9)
    ap.add_argument("--recomb", type=float, default=1e-8)
    ap.add_argument("--offset", type=int, required=True,
                    help="Generations between T1 and T2")
    ap.add_argument("--no_recent_mutations", action="store_true",
                    help="Skip overlaying mutations on the forward-simulated period")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rep", type=int, default=None)
    ap.add_argument("--sel_s", type=float, default=None)
    ap.add_argument("--decline_rate", type=float, default=None)
    ap.add_argument("--f0", type=float, default=None)
    ap.add_argument("--tsv_out", default=None)
    ap.add_argument("--out_trees", default=None)
    args = ap.parse_args()

    ts = tskit.load(args.trees)

    # Safety net: should never trigger, because the msprime ancestry was complete
    recapitated = False
    if max_roots(ts) > 1:
        ts = pyslim.recapitate(ts, ancestral_Ne=args.Ne,
                               recombination_rate=args.recomb,
                               random_seed=args.seed)
        recapitated = True

    if not args.no_recent_mutations:
        # Forward period spans tskit times 0 .. offset (T1 samples sit at time = offset)
        ts = msprime.sim_mutations(
            ts,
            rate=args.mu,
            model=SLIM_MUT_MODEL(type=0, next_id=pyslim.next_slim_mutation_id(ts)),
            start_time=0,
            end_time=args.offset,
            keep=True,
            random_seed=args.seed,
        )

    if args.out_trees:
        ts.dump(args.out_trees)

    with open(args.site_file) as f:
        sel_pos, f0_realized, attempts, init_seed = f.read().split()[:4]
    sel_pos = int(sel_pos)

    n1 = sample_nodes_from_pedigrees(ts, load_ids(args.t1_ids))
    n2 = sample_nodes_from_pedigrees(ts, load_ids(args.t2_ids))

    fst, pi1, pi2, dxy = hudson_fst(ts, n1, n2)
    D1 = tajimas_d(ts, n1)
    D2 = tajimas_d(ts, n2)
    p1 = site_freq(ts, sel_pos, n1)
    p2 = site_freq(ts, sel_pos, n2)

    print(f"nodes_t1={len(n1)} nodes_t2={len(n2)} sites={ts.num_sites} "
          f"recapitated={recapitated}")
    print(f"pi_t1={pi1:.6g} pi_t2={pi2:.6g} dxy={dxy:.6g} FST={fst:.6g}")
    print(f"TajD_t1={D1:.6g} TajD_t2={D2:.6g} delta={D2 - D1:.6g}")
    print(f"selected site {sel_pos}: popF0={f0_realized} sample p_t1={p1:.4f} p_t2={p2:.4f}")

    if args.tsv_out:
        header = ["rep", "sel_s", "decline_rate", "offset", "model", "seed",
                  "num_sites", "num_mutations",
                  "pi_t1", "pi_t2", "dxy", "fst_hudson",
                  "tajd_t1", "tajd_t2", "delta_tajd",
                  # new columns
                  "f0_target", "f0_realized", "sel_pos", "p_sel_t1", "p_sel_t2",
                  "init_attempts", "recapitated"]
        row = [args.rep, args.sel_s, args.decline_rate, args.offset, "msprime_init",
               args.seed, ts.num_sites, ts.num_mutations,
               pi1, pi2, dxy, fst, D1, D2, D2 - D1,
               args.f0, f0_realized, sel_pos, p1, p2, attempts, recapitated]
        new = not os.path.exists(args.tsv_out)
        with open(args.tsv_out, "a") as f:
            if new:
                f.write("\t".join(header) + "\n")
            f.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
