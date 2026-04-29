#!/usr/bin/env python3
"""
parse_kraken2_reports.py
------------------------
Parse a directory of Kraken2 report files and produce a single combined
abundance table suitable for downstream R analysis (e.g. phyloseq).

Usage
-----
    python parse_kraken2_reports.py \\
        --report-dir  /path/to/kraken2_reports \\
        --rank        S \\
        --min-reads   10 \\
        --output      combined_kraken2_species.tsv

Output columns: taxon_name, taxon_id, <sample_1>, <sample_2>, ...
"""

import argparse
import logging
import os
import sys
from pathlib import Path

import pandas as pd

# Kraken2 report column positions (standard format)
K2_COLS = ["pct", "clade_reads", "taxon_reads", "rank", "taxid", "name"]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def parse_k2_report(path: Path, rank: str) -> pd.Series:
    """Return a Series of read counts for the requested taxonomic rank."""
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=K2_COLS,
        dtype={"taxid": str},
    )
    df["name"] = df["name"].str.strip()

    # Select exact rank rows (e.g. "S" for species, "G" for genus)
    subset = df[df["rank"] == rank].copy()
    subset = subset.set_index("name")["taxon_reads"]
    return subset


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--report-dir",  required=True,
                        help="Directory containing *.k2report files")
    parser.add_argument("--rank",        default="S",
                        help="Taxonomic rank code: D, P, C, O, F, G, S (default: S)")
    parser.add_argument("--min-reads",   type=int, default=0,
                        help="Minimum total reads across all samples to retain a taxon")
    parser.add_argument("--output",      required=True,
                        help="Output TSV file path")
    args = parser.parse_args(argv)

    report_dir = Path(args.report_dir)
    report_files = sorted(report_dir.glob("*.k2report"))

    if not report_files:
        log.error("No *.k2report files found in: %s", report_dir)
        return 1

    log.info("Found %d report files in %s", len(report_files), report_dir)

    abundance_dict: dict[str, pd.Series] = {}
    for rf in report_files:
        sample_id = rf.stem  # filename without extension
        try:
            abundance_dict[sample_id] = parse_k2_report(rf, args.rank)
            log.info("  Parsed %-30s  taxa at rank '%s': %d",
                     rf.name, args.rank, len(abundance_dict[sample_id]))
        except Exception as exc:  # noqa: BLE001
            log.warning("  Skipping %s — error: %s", rf.name, exc)

    if not abundance_dict:
        log.error("No samples were successfully parsed.")
        return 1

    combined = (
        pd.DataFrame(abundance_dict)
        .fillna(0)
        .astype(int)
    )

    # Apply minimum-read filter
    if args.min_reads > 0:
        before = len(combined)
        combined = combined[combined.sum(axis=1) >= args.min_reads]
        log.info("Low-count filter (< %d reads): removed %d taxa, kept %d",
                 args.min_reads, before - len(combined), len(combined))

    combined.index.name = "taxon_name"
    combined = combined.sort_values(combined.columns.tolist(), ascending=False)

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, sep="\t")

    log.info("Written: %s  (%d taxa × %d samples)",
             out_path, *combined.shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
