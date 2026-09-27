#!/usr/bin/env python3
import argparse
import warnings
import json
import time
import os
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
    """Coerce scalar/length-1 tskit statistic output to a Python float."""
    arr = np.asarray(x)
    return float(arr) if arr.shape == () else float(arr.ravel()[0])


def nucleotide_diversity(ts, samples):
    """Mean pairwise nucleotide diversity per unit sequence length."""
    return _as_float(ts.diversity(sample_sets=[samples]))


def get_ts_population_name(ts, j):
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


def allele_frequency_trajectories(ts, nodes_t1, nodes_t2, nodes_t3):
    """
    Compute total derived-allele frequency at every site for t1, t2, and t3.

    With the default tskit allele encoding, genotype 0 is the ancestral state
    and values >0 are derived states. We collapse all derived states at a site
    into one total derived-allele frequency. This is appropriate here because
    recurrent/multiallelic mutation at the same position should be rare.

    Sites that are ancestral (frequency 0) in all three temporal samples are
    omitted. Such sites can exist in the full tree sequence because a mutation
    may occur only on lineages outside these remembered samples.
    """
    n1 = len(nodes_t1)
    n2 = len(nodes_t2)
    n3 = len(nodes_t3)

    all_nodes = np.concatenate([nodes_t1, nodes_t2, nodes_t3]).astype(np.int32)

    positions = []
    af1 = []
    af2 = []
    af3 = []

    for var in ts.variants(samples=all_nodes, isolated_as_missing=False):
        g = var.genotypes

        g1 = g[:n1]
        g2 = g[n1:n1 + n2]
        g3 = g[n1 + n2:n1 + n2 + n3]

        # var.alleles[0] is ancestral under tskit's default encoding.
        p1 = np.mean(g1 != 0)
        p2 = np.mean(g2 != 0)
        p3 = np.mean(g3 != 0)

        # Exclude mutations absent from all three temporal samples.
        if p1 == 0.0 and p2 == 0.0 and p3 == 0.0:
            continue

        positions.append(var.site.position)
        af1.append(p1)
        af2.append(p2)
        af3.append(p3)

    return (
        np.asarray(positions, dtype=float),
        np.asarray(af1, dtype=float),
        np.asarray(af2, dtype=float),
        np.asarray(af3, dtype=float),
    )


def summarize_delta_af(af_a, af_b):
    """
    Summarize allele-frequency change from sample A -> B.

    Returns:
      mean_delta_af:
          mean(p_B - p_A), i.e. signed mean change across sites

      abs_mean_delta_af:
          abs(mean(p_B - p_A)); useful if your downstream statistic is the
          absolute value of a window/segment-level mean delta AF

      mean_abs_delta_af:
          mean(abs(p_B - p_A)); magnitude of sitewise allele-frequency change

    These differ because positive and negative sitewise changes can cancel in
    mean_delta_af before taking the absolute value.
    """
    if af_a.size == 0:
        return np.nan, np.nan, np.nan

    delta = af_b - af_a
    mean_delta = float(np.mean(delta))
    return (
        mean_delta,
        abs(mean_delta),
        float(np.mean(np.abs(delta))),
    )


def write_site_af_tsv(path, positions, af1, af2, af3):
    """Write per-site temporal allele frequencies and changes."""
    with open(path, "w") as f:
        f.write(
            "position\taf_t1\taf_t2\taf_t3"
            "\tdelta_af_t1_t2\tdelta_af_t2_t3\tdelta_af_t1_t3\n"
        )
        for pos, p1, p2, p3 in zip(positions, af1, af2, af3):
            f.write(
                f"{pos}\t{p1}\t{p2}\t{p3}"
                f"\t{p2 - p1}\t{p3 - p2}\t{p3 - p1}\n"
            )


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Recapitate + overlay neutral mutations, then compute nucleotide "
            "diversity and temporal allele-frequency change across t1, t2, and t3."
        )
    )

    ap.add_argument(
        "--tsv_out",
        default=None,
        help="If set, append a one-line TSV of summary statistics to this file.",
    )
    ap.add_argument(
        "--site_tsv_out",
        default=None,
        help="If set, write per-site AF at t1/t2/t3 and pairwise delta AF.",
    )

    ap.add_argument("--rep", type=int, default=None, help="Replicate ID for reporting.")
    ap.add_argument("--offset2", type=int, default=None, help="Generations from t1 to t2, for reporting.")
    ap.add_argument("--offset3", type=int, default=None, help="Generations from t1 to t3, for reporting.")
    ap.add_argument("--sel_s", type=float, default=None, help="Selection coefficient for reporting.")
    ap.add_argument("--decline_rate", type=float, default=None, help="Decline rate for reporting.")

    ap.add_argument("--trees", required=True)
    ap.add_argument("--t1_ids", required=True)
    ap.add_argument("--t2_ids", required=True)
    ap.add_argument("--t3_ids", required=True)

    ap.add_argument("--Ne", type=float, default=150_000)
    ap.add_argument("--mu", type=float, default=2.04e-9)
    ap.add_argument("--recomb", type=float, default=1.0e-8)
    ap.add_argument(
        "--L",
        type=int,
        default=None,
        help=(
            "Sequence length. If omitted, use the length stored in the tree sequence. "
            "If supplied, it must match the tree-sequence length."
        ),
    )
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
    log(
        f"      Loaded in {toc(t0):.2f}s | L={ts.sequence_length} | "
        f"pops={ts.num_populations} | inds={ts.num_individuals} | "
        f"samples={ts.num_samples}",
        args.verbose,
    )

    if args.debug_metadata:
        log(
            "      DEBUG: population names (inferred): "
            + ", ".join(get_ts_population_name(ts, j) for j in range(ts.num_populations)),
            True,
        )
        if ts.num_individuals > 0:
            md0 = ts.individual(0).metadata
            log(f"      DEBUG: individual[0].metadata type={type(md0)}", True)
            md0d = _coerce_metadata_to_dict(md0)
            log(f"      DEBUG: individual[0].metadata parsed dict={md0d}", True)

    if args.L is None:
        sequence_length = ts.sequence_length
    else:
        sequence_length = float(args.L)
        if not np.isclose(sequence_length, ts.sequence_length):
            raise ValueError(
                f"--L ({args.L}) does not match tree-sequence length "
                f"({ts.sequence_length}). Use the SLiM segment length or omit --L."
            )

    # ---- Demography ----
    log("[2/7] Building demography", args.verbose)
    t0 = tic()
    demography = None
    pop_size = args.Ne
    if ts.num_populations > 1:
        demography = build_demography_from_ts(ts, args.Ne)
        pop_size = None
        log(
            f"      Using demography with names: {[p.name for p in demography.populations]}",
            args.verbose,
        )
    log(f"      Demography ready in {toc(t0):.2f}s", args.verbose)

    # ---- Recapitate ----
    log("[3/7] Recapitating (this is usually the slow step)...", args.verbose)
    t0 = tic()
    recomb_map = msprime.RateMap.uniform(
        sequence_length=sequence_length,
        rate=args.recomb,
    )

    ts_recap = msprime.sim_ancestry(
        initial_state=ts,
        recombination_rate=recomb_map,
        population_size=pop_size,
        demography=demography,
        model=args.model,
    )
    log(
        f"      Recapitation done in {toc(t0):.2f}s | "
        f"trees={ts_recap.num_trees} | nodes={ts_recap.num_nodes}",
        args.verbose,
    )

    # ---- Mutation overlay ----
    log("[4/7] Overlaying neutral mutations", args.verbose)
    t0 = tic()
    ts_mut = msprime.sim_mutations(
        ts_recap,
        rate=args.mu,
        model=msprime.SLiMMutationModel(type=0),
        keep=True,
        random_seed=args.seed,
    )
    log(
        f"      Mutation overlay done in {toc(t0):.2f}s | "
        f"sites={ts_mut.num_sites} | muts={ts_mut.num_mutations}",
        args.verbose,
    )

    if args.out_trees is not None:
        log(f"[5/7] Writing recap+mut trees to: {args.out_trees}", args.verbose)
        t0 = tic()
        ts_mut.dump(args.out_trees)
        log(f"      Wrote in {toc(t0):.2f}s", args.verbose)
    else:
        log("[5/7] Skipping write (no --out_trees)", args.verbose)

    # ---- Map temporal samples ----
    log("[6/7] Mapping pedigree IDs -> nodes", args.verbose)
    t0 = tic()
    ids_t1 = load_ids(args.t1_ids)
    ids_t2 = load_ids(args.t2_ids)
    ids_t3 = load_ids(args.t3_ids)

    nodes_t1 = sample_nodes_from_pedigrees(ts_mut, ids_t1)
    nodes_t2 = sample_nodes_from_pedigrees(ts_mut, ids_t2)
    nodes_t3 = sample_nodes_from_pedigrees(ts_mut, ids_t3)

    log(
        f"      Mapped in {toc(t0):.2f}s | "
        f"nodes_t1={len(nodes_t1)} | nodes_t2={len(nodes_t2)} | nodes_t3={len(nodes_t3)}",
        args.verbose,
    )

    # ---- Statistics ----
    log("[7/7] Computing pi, delta pi, and temporal delta AF", args.verbose)
    t0 = tic()

    # Nucleotide diversity at each sampling point
    pi1 = nucleotide_diversity(ts_mut, nodes_t1)
    pi2 = nucleotide_diversity(ts_mut, nodes_t2)
    pi3 = nucleotide_diversity(ts_mut, nodes_t3)

    delta_pi_12 = pi2 - pi1
    delta_pi_23 = pi3 - pi2
    delta_pi_13 = pi3 - pi1

    # Allele-frequency trajectories at the same sites in all three samples
    positions, af1, af2, af3 = allele_frequency_trajectories(
        ts_mut,
        nodes_t1,
        nodes_t2,
        nodes_t3,
    )

    mean_daf_12, abs_mean_daf_12, mean_abs_daf_12 = summarize_delta_af(af1, af2)
    mean_daf_23, abs_mean_daf_23, mean_abs_daf_23 = summarize_delta_af(af2, af3)
    mean_daf_13, abs_mean_daf_13, mean_abs_daf_13 = summarize_delta_af(af1, af3)

    if args.site_tsv_out is not None:
        write_site_af_tsv(args.site_tsv_out, positions, af1, af2, af3)

    log(f"      Stats computed in {toc(t0):.2f}s", args.verbose)

    # ---- Human-readable output ----
    print("=== Inputs ===")
    print("trees:", args.trees)
    print("t1_ids:", args.t1_ids, " (n pedigree IDs:", len(ids_t1), "; n sample nodes:", len(nodes_t1), ")")
    print("t2_ids:", args.t2_ids, " (n pedigree IDs:", len(ids_t2), "; n sample nodes:", len(nodes_t2), ")")
    print("t3_ids:", args.t3_ids, " (n pedigree IDs:", len(ids_t3), "; n sample nodes:", len(nodes_t3), ")")
    print("")

    print("=== Recap + overlay ===")
    print(
        "model:", args.model,
        "Ne:", args.Ne,
        "mu:", args.mu,
        "recomb:", args.recomb,
        "L:", sequence_length,
    )
    print(
        "pops:", ts_mut.num_populations,
        "sites:", ts_mut.num_sites,
        "mutations:", ts_mut.num_mutations,
    )
    if args.out_trees:
        print("wrote recap+mut trees:", args.out_trees)
    if args.site_tsv_out:
        print("wrote site-level AF file:", args.site_tsv_out)
    print("")

    print("=== Nucleotide diversity ===")
    print(f"pi(t1) = {pi1:.6g}")
    print(f"pi(t2) = {pi2:.6g}")
    print(f"pi(t3) = {pi3:.6g}")
    print(f"delta_pi(t1->t2) = {delta_pi_12:.6g}")
    print(f"delta_pi(t2->t3) = {delta_pi_23:.6g}")
    print(f"delta_pi(t1->t3) = {delta_pi_13:.6g}")
    print("")

    print("=== Allele-frequency change ===")
    print("sites represented in >=1 temporal sample:", len(positions))
    print("")
    print("Signed mean delta AF:")
    print(f"  mean_delta_AF(t1->t2) = {mean_daf_12:.6g}")
    print(f"  mean_delta_AF(t2->t3) = {mean_daf_23:.6g}")
    print(f"  mean_delta_AF(t1->t3) = {mean_daf_13:.6g}")
    print("")
    print("Absolute value of segment-level mean delta AF:")
    print(f"  abs_mean_delta_AF(t1->t2) = {abs_mean_daf_12:.6g}")
    print(f"  abs_mean_delta_AF(t2->t3) = {abs_mean_daf_23:.6g}")
    print(f"  abs_mean_delta_AF(t1->t3) = {abs_mean_daf_13:.6g}")
    print("")
    print("Mean sitewise absolute delta AF:")
    print(f"  mean_abs_delta_AF(t1->t2) = {mean_abs_daf_12:.6g}")
    print(f"  mean_abs_delta_AF(t2->t3) = {mean_abs_daf_23:.6g}")
    print(f"  mean_abs_delta_AF(t1->t3) = {mean_abs_daf_13:.6g}")
    print("")
    print("offset2:", args.offset2)
    print("offset3:", args.offset3)

    # ---- One-line summary TSV ----
    if args.tsv_out is not None:
        header = "\t".join([
            "rep", "sel_s", "decline_rate", "offset2", "offset3", "model", "seed",
            "num_sites", "num_mutations", "num_af_sites",
            "pi_t1", "pi_t2", "pi_t3",
            "delta_pi_t1_t2", "delta_pi_t2_t3", "delta_pi_t1_t3",
            "mean_delta_af_t1_t2", "mean_delta_af_t2_t3", "mean_delta_af_t1_t3",
            "abs_mean_delta_af_t1_t2", "abs_mean_delta_af_t2_t3", "abs_mean_delta_af_t1_t3",
            "mean_abs_delta_af_t1_t2", "mean_abs_delta_af_t2_t3", "mean_abs_delta_af_t1_t3",
        ])

        row = "\t".join(map(str, [
            args.rep, args.sel_s, args.decline_rate, args.offset2, args.offset3, args.model, args.seed,
            ts_mut.num_sites, ts_mut.num_mutations, len(positions),
            pi1, pi2, pi3,
            delta_pi_12, delta_pi_23, delta_pi_13,
            mean_daf_12, mean_daf_23, mean_daf_13,
            abs_mean_daf_12, abs_mean_daf_23, abs_mean_daf_13,
            mean_abs_daf_12, mean_abs_daf_23, mean_abs_daf_13,
        ]))

        if not os.path.exists(args.tsv_out):
            with open(args.tsv_out, "w") as f:
                f.write(header + "\n")

        with open(args.tsv_out, "a") as f:
            f.write(row + "\n")


if __name__ == "__main__":
    main()