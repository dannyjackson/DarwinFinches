#!/usr/bin/env python3
import argparse
import warnings
import json
import time
import numpy as np
import msprime
import tskit


def log(msg, verbose=True):
    if verbose:
        print(msg, flush=True)


def tic():
    return time.time()


def toc(t0):
    return time.time() - t0


def load_ids(path):
    with open(path) as f:
        return set(int(line.strip()) for line in f if line.strip())


def _coerce_metadata_to_dict(md):
    if md is None:
        return None
    if isinstance(md, dict):
        return md
    if isinstance(md, (bytes, bytearray, memoryview)):
        b = bytes(md)
        if len(b) == 0:
            return None
        try:
            return json.loads(b.decode("utf-8"))
        except Exception:
            return None
    if isinstance(md, str):
        s = md.strip()
        if not s:
            return None
        try:
            return json.loads(s)
        except Exception:
            return None
    return None


def get_pedigree_id(ind):
    md = _coerce_metadata_to_dict(ind.metadata)
    if not isinstance(md, dict):
        return None
    for key in ("pedigree_id", "pedigreeID", "pedigreeId"):
        if key in md:
            return md[key]
    return None


def sample_nodes_from_pedigrees(ts, pedigree_ids):
    nodes = []
    matched = 0
    for ind in ts.individuals():
        pid = get_pedigree_id(ind)
        if pid in pedigree_ids:
            matched += 1
            for n in ind.nodes:
                if n != tskit.NULL:
                    nodes.append(n)
    nodes = np.array(nodes, dtype=np.int32)
    if nodes.size == 0:
        raise RuntimeError(
            "No nodes found for the provided pedigree IDs.\n"
            f"  pedigree IDs provided: {len(pedigree_ids)}\n"
            f"  individuals matched by pedigree: {matched}\n\n"
            "Debug tip: run with --debug_metadata to print example metadata."
        )
    return nodes


def _as_float(x):
    """
    tskit stats return shapes differ by version:
      - sometimes a numpy array
      - sometimes a scalar float/np.float64
    This coerces either to a plain float.
    """
    arr = np.asarray(x)
    return float(arr) if arr.shape == () else float(arr.ravel()[0])


def hudson_fst(ts, samples_a, samples_b):
    """
    Hudson-style FST:
      FST = 1 - pi_within / dxy
      pi_within = (pi_a + pi_b)/2

    Uses:
      pi_* via ts.diversity
      dxy via ts.divergence
    """
    pi_a = _as_float(ts.diversity(sample_sets=[samples_a]))
    pi_b = _as_float(ts.diversity(sample_sets=[samples_b]))
    dxy  = _as_float(ts.divergence(sample_sets=[samples_a, samples_b]))

    pi_within = 0.5 * (pi_a + pi_b)

    if dxy == 0.0:
        return np.nan, pi_a, pi_b, dxy

    fst = 1.0 - (pi_within / dxy)
    return fst, pi_a, pi_b, dxy


def tajimas_d(ts, samples):
    if hasattr(ts, "Tajimas_D"):
        return ts.Tajimas_D(sample_sets=[samples])[0]

    n = len(samples)
    if n < 2:
        return np.nan

    S = 0
    for var in ts.variants(samples=samples, alleles="ACGT"):
        g = var.genotypes
        g = g[g >= 0]
        if g.size == 0:
            continue
        if np.any(g != g[0]):
            S += 1

    pi = ts.diversity(sample_sets=[samples])[0]

    a1 = np.sum(1.0 / np.arange(1, n))
    a2 = np.sum(1.0 / (np.arange(1, n) ** 2))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1**2))
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)

    pi_total = pi * ts.sequence_length
    thetaW_total = (S / a1)

    var = (e1 * S) + (e2 * S * (S - 1))
    if var <= 0:
        return np.nan
    return (pi_total - thetaW_total) / np.sqrt(var)


def get_ts_population_name(ts, j):
    """
    Your TS appears to use 'pop_0', 'pop_1', ... (from the earlier error).
    We try metadata first; otherwise fall back to pop_{j}.
    """
    pop = ts.population(j)
    md = _coerce_metadata_to_dict(pop.metadata)
    if isinstance(md, dict):
        for key in ("name", "pop_name", "population_name"):
            if key in md and isinstance(md[key], str) and md[key]:
                return md[key]
    return f"pop_{j}"


def build_demography_from_ts(ts, Ne):
    demography = msprime.Demography()
    for j in range(ts.num_populations):
        demography.add_population(name=get_ts_population_name(ts, j), initial_size=Ne)
    return demography


def main():
    ap = argparse.ArgumentParser(
        description="Recapitate + overlay neutral muts, then compute FST and ΔTajima's D between t1 and t2."
    )
    ap.add_argument("--tsv_out", default=None, help="If set, append a one-line TSV of summary stats to this file.")
    ap.add_argument("--rep", type=int, default=None, help="Replicate ID for reporting.")
    ap.add_argument("--sel_s", type=float, default=None, help="Selection coefficient for reporting.")
    ap.add_argument("--decline_rate", type=float, default=None, help="Decline rate for reporting.")
    ap.add_argument("--trees", required=True)
    ap.add_argument("--t1_ids", required=True)
    ap.add_argument("--t2_ids", required=True)
    ap.add_argument("--Ne", type=float, default=150_000)
    ap.add_argument("--mu", type=float, default=2.04e-9)
    ap.add_argument("--recomb", type=float, default=1.0e-8)
    ap.add_argument("--L", type=int, default=50_000)
    ap.add_argument("--model", default="dtwf", choices=["dtwf", "hudson"])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out_trees", default=None)

    ap.add_argument("--suppress_time_warning", action="store_true")
    ap.add_argument("--debug_metadata", action="store_true")
    ap.add_argument("--verbose", action="store_true", help="Print progress + timings")

    args = ap.parse_args()

    if args.suppress_time_warning:
        warnings.simplefilter("ignore", msprime.TimeUnitsMismatchWarning)

    # ---- Load ----
    log(f"[1/7] Loading tree sequence: {args.trees}", args.verbose)
    t0 = tic()
    ts = tskit.load(args.trees)
    log(f"      Loaded in {toc(t0):.2f}s | L={ts.sequence_length} | pops={ts.num_populations} | inds={ts.num_individuals} | samples={ts.num_samples}",
        args.verbose)

    if args.debug_metadata:
        log("      DEBUG: population names (inferred): " +
            ", ".join(get_ts_population_name(ts, j) for j in range(ts.num_populations)),
            True)
        if ts.num_individuals > 0:
            md0 = ts.individual(0).metadata
            log(f"      DEBUG: individual[0].metadata type={type(md0)}", True)
            md0d = _coerce_metadata_to_dict(md0)
            log(f"      DEBUG: individual[0].metadata parsed dict={md0d}", True)

    # ---- Demography ----
    log("[2/7] Building demography", args.verbose)
    t0 = tic()
    demography = None
    pop_size = args.Ne
    if ts.num_populations > 1:
        demography = build_demography_from_ts(ts, args.Ne)
        pop_size = None  # can't use population_size when demography is provided
        log(f"      Using demography with names: {[p.name for p in demography.populations]}", args.verbose)
    log(f"      Demography ready in {toc(t0):.2f}s", args.verbose)

    # ---- Recapitate ----
    log("[3/7] Recapitating (this is usually the slow step)...", args.verbose)
    t0 = tic()
    recomb_map = msprime.RateMap.uniform(sequence_length=args.L, rate=args.recomb)

    ts_recap = msprime.sim_ancestry(
        initial_state=ts,
        recombination_rate=recomb_map,
        population_size=pop_size,
        demography=demography,
        model=args.model,
    )
    log(f"      Recapitation done in {toc(t0):.2f}s | trees={ts_recap.num_trees} | nodes={ts_recap.num_nodes}",
        args.verbose)

    # ---- Mut overlay ----
    log("[4/7] Overlaying neutral mutations", args.verbose)
    t0 = tic()
    ts_mut = msprime.sim_mutations(
        ts_recap,
        rate=args.mu,
        model=msprime.SLiMMutationModel(type=0),
        keep=True,
        random_seed=args.seed,
    )
    log(f"      Mut overlay done in {toc(t0):.2f}s | sites={ts_mut.num_sites} | muts={ts_mut.num_mutations}",
        args.verbose)

    if args.out_trees is not None:
        log(f"[5/7] Writing recap+mut trees to: {args.out_trees}", args.verbose)
        t0 = tic()
        ts_mut.dump(args.out_trees)
        log(f"      Wrote in {toc(t0):.2f}s", args.verbose)
    else:
        log("[5/7] Skipping write (no --out_trees)", args.verbose)

    # ---- Map samples ----
    log("[6/7] Mapping pedigree IDs -> nodes", args.verbose)
    t0 = tic()
    ids_t1 = load_ids(args.t1_ids)
    ids_t2 = load_ids(args.t2_ids)
    nodes_t1 = sample_nodes_from_pedigrees(ts_mut, ids_t1)
    nodes_t2 = sample_nodes_from_pedigrees(ts_mut, ids_t2)
    log(f"      Mapped in {toc(t0):.2f}s | nodes_t1={len(nodes_t1)} | nodes_t2={len(nodes_t2)}",
        args.verbose)

    # ---- Stats ----
    log("[7/7] Computing statistics (FST, Tajima's D, ΔD)", args.verbose)
    t0 = tic()
    ts_t1 = ts_mut.simplify(samples=nodes_t1, keep_unary=True)
    ts_t2 = ts_mut.simplify(samples=nodes_t2, keep_unary=True)

    D_t1 = tajimas_d(ts_t1, np.arange(ts_t1.num_samples, dtype=np.int32))
    D_t2 = tajimas_d(ts_t2, np.arange(ts_t2.num_samples, dtype=np.int32))
    delta_D = D_t2 - D_t1

    fst, pi1, pi2, dxy = hudson_fst(ts_mut, nodes_t1, nodes_t2)
    log(f"      Stats computed in {toc(t0):.2f}s", args.verbose)

    # ---- Output ----
    print("=== Inputs ===")
    print("trees:", args.trees)
    print("t1_ids:", args.t1_ids, " (n pedigree IDs:", len(ids_t1), "; n sample nodes:", len(nodes_t1), ")")
    print("t2_ids:", args.t2_ids, " (n pedigree IDs:", len(ids_t2), "; n sample nodes:", len(nodes_t2), ")")
    print("")
    print("=== Recap + overlay ===")
    print("model:", args.model, "Ne:", args.Ne, "mu:", args.mu, "recomb:", args.recomb, "L:", args.L)
    print("pops:", ts_mut.num_populations, "sites:", ts_mut.num_sites, "mutations:", ts_mut.num_mutations)
    if args.out_trees:
        print("wrote recap+mut trees:", args.out_trees)
    print("")
    print("=== Results ===")
    print(f"pi(t1) = {pi1:.6g}")
    print(f"pi(t2) = {pi2:.6g}")
    print(f"dxy(t1,t2) = {dxy:.6g}")
    print(f"FST_Hudson(t1,t2) = {fst:.6g}")
    print(f"TajimaD(t1) = {D_t1:.6g}")
    print(f"TajimaD(t2) = {D_t2:.6g}")
    print(f"delta_TajimaD (t2 - t1) = {delta_D:.6g}")

    if args.tsv_out is not None:
        header = "\t".join([
            "rep", "sel_s", "decline_rate", "model", "seed",
            "num_sites", "num_mutations",
            "pi_t1", "pi_t2", "dxy", "fst_hudson",
            "tajd_t1", "tajd_t2", "delta_tajd"
        ])
        row = "\t".join(map(str, [
            args.rep, args.sel_s, args.decline_rate, args.model, args.seed,
            ts_mut.num_sites, ts_mut.num_mutations,
            pi1, pi2, dxy, fst,
            D_t1, D_t2, delta_D
        ]))
        # write header if file doesn't exist yet
        import os
        if not os.path.exists(args.tsv_out):
            with open(args.tsv_out, "w") as f:
                f.write(header + "\n")
        with open(args.tsv_out, "a") as f:
            f.write(row + "\n")


if __name__ == "__main__":
    main()
