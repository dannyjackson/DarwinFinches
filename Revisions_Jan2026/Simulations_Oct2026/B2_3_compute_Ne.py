#!/usr/bin/env python3
"""
B2_3_compute_Ne.py

Effective population size from Watterson's theta (ANGSD thetaStat do_stat
.pestPG files), with the mutation rate converted to a PER-GENERATION rate.

    theta_W (per site) = sum(tW) / sum(nSites)        over autosomes
    mu_gen             = mu_year * generation_time
    Ne                 = theta_W / (4 * mu_gen)

Why the conversion: the Darwin's finch rate of 2.04e-9 (Lamichhaney et al.
2015) is a substitution rate PER SITE PER YEAR (synonymous divergence since
the split from zebra finch). theta = 4 Ne mu needs mu per generation, so the
yearly rate is multiplied by the generation time. The previous version
divided by 4 * 2.04e-9 directly, which returns Ne x generation time.

Other fixes relative to the previous B2 notes:
  * every population uses its own sum(nSites) as the denominator
    (previously crapre used sum(nSites) and the other five a hard-coded
    946,505,819);
  * all six populations are processed (crapost was dropped from the loop);
  * species values are the mean of the pre and post estimates, computed
    from the per-population numbers (the old "Results" block mislabelled
    pre/post means as post values);
  * standard errors from a weighted delete-one-chromosome jackknife
    (Busing et al. 1999), since chromosomes differ greatly in length.

Uses only the Python 3 standard library (no environment needed).

Usage (from the effective_pop directory):
    python3 B2_3_compute_Ne.py \
        --theta_dir thetas \
        --out Ne_estimates.tsv

Options:
    --pops        populations (default: crapre crapost forpre forpost parpre parpost)
    --mu_year     per-site per-year mutation rate (default 2.04e-9)
    --gen_time    primary generation time in years (default 5)
    --gen_sens    extra generation times for a sensitivity table
                  (default 3 4.5 7; 4.5 is the value used to scale MSMC)
    --exclude     chromosomes to drop (default NC_044601, the Z)
    --prefix      keep only chromosomes starting with this (default NC_,
                  i.e. drop unplaced scaffolds)
"""
import argparse
import math
import os
import sys


def read_pestpg(path, prefix, exclude):
    """Return list of (chrom, sum_tW, nSites) for kept chromosomes."""
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        # header's first field is '#(indexStart,indexStop)(...)(WinStart,WinStop)'
        cols = {name.lstrip("#"): i for i, name in enumerate(header)}
        for need in ("Chr", "tW", "nSites"):
            if need not in cols:
                sys.exit(f"ERROR: column '{need}' not found in {path}; header = {header}")
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            chrom = p[cols["Chr"]]
            if prefix and not chrom.startswith(prefix):
                continue
            # match with or without the version suffix (NC_044601 drops NC_044601.1)
            if any(chrom == e or chrom.startswith(e + ".") for e in exclude):
                continue
            tw = float(p[cols["tW"]])
            ns = float(p[cols["nSites"]])
            if ns <= 0 or not math.isfinite(tw):
                continue
            rows.append((chrom, tw, ns))
    if not rows:
        sys.exit(f"ERROR: no chromosomes kept from {path}")
    if len({r[0] for r in rows}) != len(rows):
        sys.exit(f"ERROR: duplicate chromosome names in {path}; was do_stat run with -win? "
                 "This script expects one row per chromosome.")
    return rows


def theta_and_jackknife(rows):
    """theta_W per site and weighted delete-one-block jackknife SE (Busing et al. 1999)."""
    tw_tot = sum(r[1] for r in rows)
    n_tot = sum(r[2] for r in rows)
    theta = tw_tot / n_tot
    g = len(rows)
    if g < 2:
        return theta, float("nan")
    pseudo, loo, h = [], [], []
    for _, tw, ns in rows:
        t_minus = (tw_tot - tw) / (n_tot - ns)
        hj = n_tot / ns
        loo.append(t_minus)
        h.append(hj)
        pseudo.append(hj * theta - (hj - 1.0) * t_minus)
    theta_j = g * theta - sum((1.0 - rows[j][2] / n_tot) * loo[j] for j in range(g))
    var = sum((pseudo[j] - theta_j) ** 2 / (h[j] - 1.0) for j in range(g)) / g
    return theta, math.sqrt(var)


def fmt(x, digits=0):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    return f"{x:.{digits}f}"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--theta_dir", default="thetas")
    ap.add_argument("--pops", nargs="+",
                    default=["crapre", "crapost", "forpre", "forpost", "parpre", "parpost"])
    ap.add_argument("--mu_year", type=float, default=2.04e-9)
    ap.add_argument("--gen_time", type=float, default=5.0)
    ap.add_argument("--gen_sens", type=float, nargs="*", default=[3.0, 4.5, 7.0])
    ap.add_argument("--exclude", nargs="*", default=["NC_044601"])
    ap.add_argument("--prefix", default="NC_")
    ap.add_argument("--out", default="Ne_estimates.tsv")
    args = ap.parse_args()

    gens = [args.gen_time] + [g for g in args.gen_sens if g != args.gen_time]
    exclude = set(args.exclude)

    pop_rows = []
    for pop in args.pops:
        path = os.path.join(args.theta_dir, f"{pop}.thetas.idx.pestPG")
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")
        rows = read_pestpg(path, args.prefix, exclude)
        theta, se = theta_and_jackknife(rows)
        rec = {
            "pop": pop,
            "species": pop[:3],
            "timepoint": pop[3:],
            "n_chrom": len(rows),
            "n_sites": int(sum(r[2] for r in rows)),
            "theta_W_per_site": theta,
            "theta_W_SE": se,
        }
        for g in gens:
            denom = 4.0 * args.mu_year * g
            rec[f"Ne_g{g:g}"] = theta / denom
            rec[f"Ne_g{g:g}_SE"] = se / denom
        # what the previous version reported (per-year rate used as per-generation)
        rec["Ne_legacy_mu_per_year"] = theta / (4.0 * args.mu_year)
        pop_rows.append(rec)

    # species means of pre and post (SE of a mean of two independent estimates)
    sp_rows = []
    for sp in sorted({r["species"] for r in pop_rows}):
        rs = [r for r in pop_rows if r["species"] == sp]
        rec = {"pop": f"{sp}_mean", "species": sp, "timepoint": "mean",
               "n_chrom": None, "n_sites": None,
               "theta_W_per_site": sum(r["theta_W_per_site"] for r in rs) / len(rs),
               "theta_W_SE": math.sqrt(sum(r["theta_W_SE"] ** 2 for r in rs)) / len(rs)}
        for g in gens:
            k = f"Ne_g{g:g}"
            rec[k] = sum(r[k] for r in rs) / len(rs)
            rec[k + "_SE"] = math.sqrt(sum(r[k + "_SE"] ** 2 for r in rs)) / len(rs)
        rec["Ne_legacy_mu_per_year"] = sum(r["Ne_legacy_mu_per_year"] for r in rs) / len(rs)
        sp_rows.append(rec)

    cols = ["pop", "species", "timepoint", "n_chrom", "n_sites",
            "theta_W_per_site", "theta_W_SE"]
    for g in gens:
        cols += [f"Ne_g{g:g}", f"Ne_g{g:g}_SE"]
    cols += ["Ne_legacy_mu_per_year"]

    def cell(c, v):
        if v is None:
            return "NA"
        if c.startswith("theta"):
            return f"{v:.6g}"
        if isinstance(v, float):
            return fmt(v, 0)
        return str(v)

    with open(args.out, "w") as f:
        f.write(f"# mu_year={args.mu_year:g}; Ne = theta_W / (4 * mu_year * g); "
                f"primary g = {args.gen_time:g} yr; excluded = {','.join(sorted(exclude))}; "
                f"chromosome prefix = {args.prefix}\n")
        f.write("\t".join(cols) + "\n")
        for r in pop_rows + sp_rows:
            f.write("\t".join(cell(c, r[c]) for c in cols) + "\n")

    # short screen summary
    g0 = f"Ne_g{args.gen_time:g}"
    print(f"mu = {args.mu_year:g} per site per year x {args.gen_time:g} yr = "
          f"{args.mu_year * args.gen_time:.3g} per site per generation\n")
    print(f"{'pop':<12}{'theta_W/site':>14}{'Ne (g=' + format(args.gen_time, 'g') + ')':>14}"
          f"{'± SE':>10}{'legacy Ne':>12}")
    for r in pop_rows + sp_rows:
        print(f"{r['pop']:<12}{r['theta_W_per_site']:>14.6g}{fmt(r[g0]):>14}"
              f"{fmt(r[g0 + '_SE']):>10}{fmt(r['Ne_legacy_mu_per_year']):>12}")
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
