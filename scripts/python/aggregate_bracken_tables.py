#!/usr/bin/env python3
"""
aggregate_bracken_tables.py
---------------------------
Merge per-sample Bracken output tables into a single wide-format abundance
matrix (taxa × samples) suitable for input to phyloseq, vegan, or DESeq2.

Optionally normalises to relative abundance (fractions summing to 1 per sample).

Usage
-----
    # Raw read counts
    python aggregate_bracken_tables.py \\
        --bracken-dir data/bracken/species \\
        --pattern    '*_bracken_species.txt' \\
        --output     data/bracken_species_matrix.tsv

    # Relative abundance
    python aggregate_bracken_tables.py \\
        --bracken-dir data/bracken/species \\
        --output      data/bracken_species_relab.tsv \\
        --relative
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

# Bracken output columns we care about
BRACKEN_NAME_COL   = "name"
BRACKEN_READS_COL  = "new_est_reads"
BRACKEN_RELAB_COL  = "fraction_total_reads"


def read_bracken_table(path: Path, use_relative: bool) -> pd.Series:
    """Read one Bracken output file and return a named Series."""
    df = pd.read_csv(path, sep="\t", dtype={"taxonomy_id": str})
    df[BRACKEN_NAME_COL] = df[BRACKEN_NAME_COL].str.strip()
    value_col = BRACKEN_RELAB_COL if use_relative else BRACKEN_READS_COL
    if value_col not in df.columns:
        raise KeyError(f"Column '{value_col}' not found in {path.name}. "
                       f"Available: {list(df.columns)}")
    series = df.set_index(BRACKEN_NAME_COL)[value_col]
    series.name = path.stem  # use filename stem as sample ID
    return series


def strip_suffix(series_name: str, pattern: str) -> str:
    """Remove a glob-derived suffix from a series name for cleaner column headers."""
    # e.g. 'SAMPLE_bracken_species' -> 'SAMPLE'
    for suffix in ("_bracken_species", "_bracken_genus",
                   "_bracken_family",  "_bracken_order"):
        if series_name.endswith(suffix):
            return series_name[: -len(suffix)]
    return series_name


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bracken-dir", required=True,
                        help="Directory containing Bracken output files")
    parser.add_argument("--pattern",     default="*.txt",
                        help="Glob pattern for Bracken files (default: '*.txt')")
    parser.add_argument("--output",      required=True,
                        help="Output TSV path")
    parser.add_argument("--relative",    action="store_true",
                        help="Use relative abundances instead of raw read counts")
    parser.add_argument("--min-samples", type=int, default=1,
                        help="Keep taxa present (non-zero) in at least N samples (default: 1)")
    args = parser.parse_args(argv)

    bracken_dir = Path(args.bracken_dir)
    files = sorted(bracken_dir.glob(args.pattern))

    if not files:
        log.error("No files matching '%s' found in: %s", args.pattern, bracken_dir)
        return 1

    log.info("Found %d Bracken files", len(files))

    series_list: list[pd.Series] = []
    for f in files:
        try:
            s = read_bracken_table(f, args.relative)
            s.name = strip_suffix(s.name, args.pattern)
            series_list.append(s)
            log.info("  %-35s  taxa: %d", f.name, len(s))
        except Exception as exc:  # noqa: BLE001
            log.warning("  Skipping %s — %s", f.name, exc)

    if not series_list:
        log.error("No files were successfully read.")
        return 1

    combined = pd.DataFrame(series_list).T.fillna(0)
    if not args.relative:
        combined = combined.astype(int)

    # Filter by minimum sample prevalence
    if args.min_samples > 1:
        prevalence = (combined > 0).sum(axis=1)
        before = len(combined)
        combined = combined[prevalence >= args.min_samples]
        log.info("Prevalence filter (present in ≥ %d samples): removed %d taxa, kept %d",
                 args.min_samples, before - len(combined), len(combined))

    combined.index.name = "taxon"

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, sep="\t")

    log.info("Written: %s  (%d taxa × %d samples)", out_path, *combined.shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
