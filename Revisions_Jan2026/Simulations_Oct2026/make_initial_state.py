#!/usr/bin/env python3
"""
Build an ancestry-aware starting population for SLiM.

1. Simulate the full neutral genealogy of a population of Ne diploids with
   msprime (DTWF for the recent past, then the standard coalescent).
2. Overlay neutral mutations.
3. Find a biallelic site whose population frequency is within `tol` of f0.
   If none exists, rerun steps 1-2 with a new seed (rejection sampling, so the
   genealogy is correctly conditioned on carrying a variant at ~f0).
4. Annotate the tree sequence with SLiM metadata (pyslim) so SLiM can load it
   with sim.readFromPopulationFile().
5. Write the chosen site's position and realized frequency for SLiM to read.

The chosen variant's carriers form a clade in the local gene tree, so the
allele starts in realistic linkage disequilibrium with its neighbours.

Requires: msprime >= 1.3, tskit, pyslim with SLiM 5 support (>= 1.1).
"""
import argparse
import sys

import msprime
import numpy as np
import pyslim


def simulate(Ne, L, mu, recomb, seed, dtwf_gens):
    model = [
        msprime.DiscreteTimeWrightFisher(duration=dtwf_gens),
        msprime.StandardCoalescent(),
    ]
    ts = msprime.sim_ancestry(
        samples=Ne,                 # every individual in the population
        ploidy=2,
        population_size=Ne,
        sequence_length=L,
        recombination_rate=recomb,
        model=model,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts,
        rate=mu,
        model=msprime.BinaryMutationModel(),
        discrete_genome=True,
        random_seed=seed,
    )
    return ts


def pick_site(ts, f0, tol, rng, center_frac):
    """Return (position, frequency) of a biallelic site with |freq - f0| <= tol,
    preferring the closest; ties broken at random. None if no candidate."""
    L = ts.sequence_length
    lo = L * (1 - center_frac) / 2
    hi = L - lo
    n = ts.num_samples
    cands = []
    for site in ts.sites():
        if len(site.mutations) != 1:          # skip recurrent / multiallelic sites
            continue
        if not (lo <= site.position < hi):
            continue
        m = site.mutations[0]
        if m.parent != -1:
            continue
        # number of samples below the mutation's node = derived allele count
        tree = ts.at(site.position)
        p = tree.num_samples(m.node) / n
        if abs(p - f0) <= tol:
            cands.append((abs(p - f0), site.position, p))
    if not cands:
        return None
    best = min(c[0] for c in cands)
    ties = [c for c in cands if c[0] == best]
    _, pos, p = ties[rng.integers(len(ties))]
    return int(pos), p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--Ne", type=int, required=True)
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--f0", type=float, required=True)
    ap.add_argument("--tol", type=float, default=None,
                    help="Accept sites with |freq - f0| <= tol "
                         "(default: min(0.25*f0, 0.05))")
    ap.add_argument("--mu", type=float, default=2.04e-9)
    ap.add_argument("--recomb", type=float, default=1e-8)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--dtwf_gens", type=int, default=100,
                    help="Generations simulated with DTWF before switching to "
                         "the coalescent (needed because sample size = 2Ne)")
    ap.add_argument("--center_frac", type=float, default=1.0,
                    help="Only choose sites in the central fraction of the "
                         "segment (1.0 = anywhere)")
    ap.add_argument("--max_attempts", type=int, default=50)
    ap.add_argument("--out_trees", required=True)
    ap.add_argument("--out_site", required=True,
                    help="Writes 'position<TAB>freq<TAB>attempts<TAB>seed_used'")
    args = ap.parse_args()

    tol = args.tol if args.tol is not None else min(0.25 * args.f0, 0.05)
    rng = np.random.default_rng(args.seed)

    for attempt in range(1, args.max_attempts + 1):
        seed = args.seed + 1_000_003 * (attempt - 1)   # distinct seed per attempt
        ts = simulate(args.Ne, args.L, args.mu, args.recomb, seed, args.dtwf_gens)
        hit = pick_site(ts, args.f0, tol, rng, args.center_frac)
        print(f"attempt {attempt}: seed={seed} trees={ts.num_trees} "
              f"sites={ts.num_sites} hit={hit}", flush=True)
        if hit is not None:
            break
    else:
        sys.exit(f"ERROR: no site within {tol} of f0={args.f0} after "
                 f"{args.max_attempts} attempts; increase L or --tol")

    pos, p = hit

    # SLiM metadata. tick=1, stage="late": SLiM loads this in a `1 late()` block.
    # annotate_mutations=True turns every msprime mutation into a SLiM
    # mutation of type m0 with s = 0 (the SLiM script must define m0).
    ts = pyslim.annotate(ts, model_type="WF", tick=1, stage="late")
    ts.dump(args.out_trees)

    with open(args.out_site, "w") as f:
        f.write(f"{pos}\t{p:.6f}\t{attempt}\t{seed}\n")
    print(f"chosen site pos={pos} freq={p:.4f} (target {args.f0} ± {tol})")


if __name__ == "__main__":
    main()
