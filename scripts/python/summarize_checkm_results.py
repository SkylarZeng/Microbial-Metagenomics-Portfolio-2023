#!/usr/bin/env python3
"""
summarize_checkm_results.py
---------------------------
Collect per-sample CheckM2 quality_report.tsv files from a run directory and
produce a single aggregated summary table with quality tier annotations.

Quality tiers (MIMAG standards):
  High quality    : completeness ≥ 90%  AND contamination ≤ 5%
  Medium quality  : completeness ≥ 70%  AND contamination ≤ 10%
  Low quality     : all others that passed binning

Usage
-----
    python summarize_checkm_results.py \\
        --checkm-dir /path/to/checkm2_results \\
        --output     data/checkm2_all_samples_summary.csv
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

HQ_COMPLETENESS  = 90.0
HQ_CONTAMINATION =  5.0
MQ_COMPLETENESS  = 70.0
MQ_CONTAMINATION = 10.0


def assign_quality_tier(row: pd.Series) -> str:
    """Assign MIMAG quality tier."""
    comp  = row["Completeness"]
    cont  = row["Contamination"]
    if comp >= HQ_COMPLETENESS and cont <= HQ_CONTAMINATION:
        return "high"
    if comp >= MQ_COMPLETENESS and cont <= MQ_CONTAMINATION:
        return "medium"
    return "low"


def load_quality_report(path: Path, sample_id: str) -> pd.DataFrame:
    """Load a single CheckM2 quality_report.tsv and add sample column."""
    df = pd.read_csv(path, sep="\t")
    df.insert(0, "sample_id", sample_id)
    df.rename(columns={"Name": "mag_id"}, inplace=True)
    return df


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--checkm-dir", required=True,
                        help="Root directory containing per-sample CheckM2 output folders")
    parser.add_argument("--output",     required=True,
                        help="Output CSV path for aggregated summary")
    args = parser.parse_args(argv)

    checkm_dir = Path(args.checkm_dir)
    report_files = sorted(checkm_dir.glob("*/quality_report.tsv"))

    if not report_files:
        log.error("No quality_report.tsv files found under: %s", checkm_dir)
        return 1

    log.info("Collecting CheckM2 reports from %d samples", len(report_files))

    frames: list[pd.DataFrame] = []
    for rpt in report_files:
        sample_id = rpt.parent.name
        try:
            df = load_quality_report(rpt, sample_id)
            frames.append(df)
            log.info("  %-20s  bins: %d", sample_id, len(df))
        except Exception as exc:  # noqa: BLE001
            log.warning("  Skipping %s — %s", rpt, exc)

    if not frames:
        log.error("No data loaded. Check directory structure.")
        return 1

    combined = pd.concat(frames, ignore_index=True)
    combined["quality_tier"] = combined.apply(assign_quality_tier, axis=1)

    # Summary statistics
    tier_counts = combined.groupby(["sample_id", "quality_tier"]).size().unstack(fill_value=0)
    log.info("\nQuality tier summary per sample:\n%s", tier_counts.to_string())

    overall = combined["quality_tier"].value_counts()
    log.info("\nOverall MAG counts: %s", overall.to_dict())

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_path, index=False)
    log.info("Written: %s  (%d MAGs total)", out_path, len(combined))

    # Print high/medium summary to stdout
    passing = combined[combined["quality_tier"].isin(["high", "medium"])]
    print(f"\nPassing MAGs (≥70% completeness, ≤10% contamination): {len(passing)}")
    print(f"  High quality  (≥90% / ≤5%):   {(combined['quality_tier'] == 'high').sum()}")
    print(f"  Medium quality (≥70% / ≤10%):  {(combined['quality_tier'] == 'medium').sum()}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
