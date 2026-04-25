#!/usr/bin/env python3
"""
calc_ratio.py - Calculate allele frequency (SNP-index) for each site and
                filter candidate variants using delta-SNP-index.

BSA-Seq strategy: compare mutant pool (mut) vs wild-type pool (wt).
  delta-SNP-index = mut_ratio - wt_ratio
  Sites where |delta-SNP-index| exceeds the threshold are retained as candidates.
"""

import argparse
import gzip
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser(
    description="Calculate BSA-Seq SNP-index and filter candidate variants"
)
parser.add_argument("--vcf",    required=True, help="Input VCF file (can be .gz compressed)")
parser.add_argument("--chr",    required=True, help="Target chromosome (e.g. 1 or Chr1)")
parser.add_argument("--start",  type=int, required=True, help="Region start position")
parser.add_argument("--end",    type=int, required=True, help="Region end position")
parser.add_argument("--delta",  type=float, default=0.3,
                    help="Delta-SNP-index threshold (|mut_ratio - wt_ratio| > this value), default 0.3")
parser.add_argument("--direction", choices=["both", "positive", "negative"], default="both",
                    help=("Candidate direction: both keeps |delta| sites, positive keeps "
                          "mut_ratio - wt_ratio > threshold, negative keeps the opposite. "
                          "Default: both"))
parser.add_argument("--min_dp", type=int, default=10,
                    help="Minimum sequencing depth filter (must be >= 1), default 10")
parser.add_argument("--outdir", default="./result", help="Output directory")
args = parser.parse_args()

if args.min_dp < 1:
    print("[ERROR] --min_dp must be at least 1 to prevent division by zero.")
    sys.exit(1)

os.makedirs(args.outdir, exist_ok=True)


def normalize_chrom(chrom):
    """Normalize chromosome name to support Chr1 / chr1 / 1 formats."""
    chrom = str(chrom).strip()
    if chrom.lower().startswith('chr'):
        chrom = chrom[3:]
    return chrom


target_chrom = normalize_chrom(args.chr)

records = []
opener = gzip.open if args.vcf.endswith(".gz") else open

with opener(args.vcf, "rt") as f:
    samples = []
    for line in f:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            headers = line.strip().split("\t")
            samples = headers[9:]  # Expected order: mut, wt
            continue

        cols = line.strip().split("\t")
        if len(cols) < 10:
            continue

        chrom, pos, _, ref, alt = cols[0], int(cols[1]), cols[2], cols[3], cols[4]

        if ',' in alt:
            continue

        filter_val = cols[6] if len(cols) > 6 else 'PASS'
        if filter_val not in ('PASS', '.'):
            continue

        fmt = cols[8].split(":")
        if "AD" not in fmt:
            continue
        ad_idx = fmt.index("AD")

        row = {"CHROM": chrom, "POS": pos, "REF": ref, "ALT": alt}
        valid = True
        for i, sname in enumerate(samples):
            sdata = cols[9 + i].split(":")
            if sdata[0] in ("./.", "."):
                valid = False
                break
            ad_str = sdata[ad_idx] if ad_idx < len(sdata) else "."
            if ad_str == "." or "." in ad_str.split(","):
                valid = False
                break
            ad = list(map(int, ad_str.split(",")))
            total = sum(ad)
            if total == 0 or total < args.min_dp:
                valid = False
                break
            row[f"{sname}_ratio"] = round(ad[1] / total, 4)
            row[f"{sname}_dp"]    = total

        if valid:
            records.append(row)

df = pd.DataFrame(records)
output_columns = [
    "CHROM", "POS", "REF", "ALT",
    "mut_ratio", "mut_dp", "wt_ratio", "wt_dp", "delta_index",
]

if df.empty:
    print("WARNING: No valid variant sites found. "
          "Check your VCF file, filter settings, and --min_dp threshold.")
    empty_df = pd.DataFrame(columns=output_columns)
    empty_df.to_csv(f"{args.outdir}/all_variants_ratio.csv", index=False)
    empty_df.to_csv(f"{args.outdir}/candidate_mutations.csv", index=False)
    empty_df.to_csv(f"{args.outdir}/candidate_positive_mutations.csv", index=False)
    empty_df.to_csv(f"{args.outdir}/candidate_negative_mutations.csv", index=False)
    sys.exit(0)

if "mut_ratio" not in df.columns or "wt_ratio" not in df.columns:
    print("WARNING: 'mut_ratio' or 'wt_ratio' column not found. "
          "Check that VCF sample order is: mut pool first, wt pool second.")
    df.to_csv(f"{args.outdir}/all_variants_ratio.csv", index=False)
    pd.DataFrame(columns=output_columns).to_csv(f"{args.outdir}/candidate_mutations.csv", index=False)
    pd.DataFrame(columns=output_columns).to_csv(f"{args.outdir}/candidate_positive_mutations.csv", index=False)
    pd.DataFrame(columns=output_columns).to_csv(f"{args.outdir}/candidate_negative_mutations.csv", index=False)
else:
    df["delta_index"] = (df["mut_ratio"] - df["wt_ratio"]).round(4)
    df.to_csv(f"{args.outdir}/all_variants_ratio.csv", index=False)
    print(f"Total valid variant sites: {len(df)}")

    region_mask = (
        (df["CHROM"].apply(normalize_chrom) == target_chrom) &
        (df["POS"] >= args.start) &
        (df["POS"] <= args.end)
    )

    positive = df[region_mask & (df["delta_index"] > args.delta)].sort_values(
        "delta_index", ascending=False
    )
    negative = df[region_mask & (df["delta_index"] < -args.delta)].sort_values(
        "delta_index", ascending=True
    )
    positive.to_csv(f"{args.outdir}/candidate_positive_mutations.csv", index=False)
    negative.to_csv(f"{args.outdir}/candidate_negative_mutations.csv", index=False)

    if args.direction == "positive":
        candidates = positive
    elif args.direction == "negative":
        candidates = negative
    else:
        candidates = df[
            region_mask & (df["delta_index"].abs() > args.delta)
        ].sort_values("delta_index", key=lambda x: x.abs(), ascending=False)

    candidates.to_csv(f"{args.outdir}/candidate_mutations.csv", index=False)
    print(
        f"Candidate variant sites "
        f"({args.chr}:{args.start}-{args.end}, direction={args.direction}, "
        f"threshold={args.delta}): "
        f"{len(candidates)}"
    )
    print(f"  positive candidates: {len(positive)}")
    print(f"  negative candidates: {len(negative)}")
