#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np


parser = argparse.ArgumentParser(
    description="""
    Convert a windowed genomic statistic to a standard four-column BED file:

    chromosome   start   end   score
    """
)

parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)

parser.add_argument("--chrom-col", required=True)
parser.add_argument("--score-col", required=True)

# Either midpoint...
parser.add_argument("--position-col", default=None)
parser.add_argument("--window-size", type=int, default=None)

# ...or explicit coordinates
parser.add_argument("--start-col", default=None)
parser.add_argument("--end-col", default=None)

args = parser.parse_args()


df = pd.read_csv(
    args.input,
    sep=r"\s+",
    comment="#"
)


# ------------------------------------------------------------
# Check columns
# ------------------------------------------------------------

required = [args.chrom_col, args.score_col]

for column in required:
    if column not in df.columns:
        raise ValueError(
            f"Column '{column}' not found.\n"
            f"Available columns: {list(df.columns)}"
        )


# ------------------------------------------------------------
# Get coordinates
# ------------------------------------------------------------

if args.start_col is not None and args.end_col is not None:

    if args.start_col not in df.columns:
        raise ValueError(f"Column '{args.start_col}' not found.")

    if args.end_col not in df.columns:
        raise ValueError(f"Column '{args.end_col}' not found.")

    start = pd.to_numeric(
        df[args.start_col],
        errors="coerce"
    )

    end = pd.to_numeric(
        df[args.end_col],
        errors="coerce"
    )

else:

    if args.position_col is None or args.window_size is None:
        raise ValueError(
            "Supply either:\n"
            "  --start-col and --end-col\n"
            "or:\n"
            "  --position-col and --window-size"
        )

    if args.position_col not in df.columns:
        raise ValueError(
            f"Column '{args.position_col}' not found."
        )

    midpoint = pd.to_numeric(
        df[args.position_col],
        errors="coerce"
    )

    # Convert window midpoint to BED coordinates.
    start = np.floor(
        midpoint - args.window_size / 2
    )

    start = np.maximum(start, 0)

    end = start + args.window_size


# ------------------------------------------------------------
# Extract statistic
# ------------------------------------------------------------

score = pd.to_numeric(
    df[args.score_col],
    errors="coerce"
)


out = pd.DataFrame({
    "chrom": df[args.chrom_col],
    "start": start,
    "end": end,
    "score": score
})


# Remove invalid rows
out = out.dropna(
    subset=["chrom", "start", "end", "score"]
)


out["start"] = out["start"].astype(int)
out["end"] = out["end"].astype(int)


# BED requires end > start
out = out[out["end"] > out["start"]]


# Sort
out = out.sort_values(
    ["chrom", "start", "end"]
)


# Write without header
out.to_csv(
    args.output,
    sep="\t",
    index=False,
    header=False
)


print(f"Wrote {len(out)} windows to:")
print(args.output)
