#!/usr/bin/env python3

"""
Convert a CheckM `bin_stats_ext.tsv` file into a flatter tab-delimited summary.

Usage:
    python3 checkm_summary.py bin_stats_ext.tsv bin_stats_ext.txt
"""

import json
import sys


HEADER = [
    "Bin Id",
    "Marker Lineage",
    "Genomes",
    "Markers",
    "Marker Sets",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+",
    "Completeness",
    "Contamination",
    "GC",
    "GC std",
    "Genome Size",
    "Ambiguous Bases",
    "Scaffolds",
    "Contigs",
    "Translation Table",
    "Predicted Genes",
]


def load_stats(path):
    stats = {}
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.replace("'", '"').rstrip("\n")
            bin_id, payload = line.split("\t", 1)
            normalized_id = "bin_" + bin_id.replace(".bin.", "_")
            stats[normalized_id] = json.loads(payload)
    return stats


def write_summary(stats, path):
    with open(path, "w", encoding="utf-8") as output:
        output.write("\t".join(HEADER) + "\n")
        for key, values in stats.items():
            row = [
                key,
                values["marker lineage"],
                str(values["# genomes"]),
                str(values["# markers"]),
                str(values["# marker sets"]),
                str(values["0"]),
                str(values["1"]),
                str(values["2"]),
                str(values["3"]),
                str(values["4"]),
                str(values["5+"]),
                str(values["Completeness"]),
                str(values["Contamination"]),
                str(values["GC"]),
                str(values["GC std"]),
                str(values["Genome size"]),
                str(values["# ambiguous bases"]),
                str(values["# scaffolds"]),
                str(values["# contigs"]),
                str(values["Translation table"]),
                str(values["# predicted genes"]),
            ]
            output.write("\t".join(row) + "\n")


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python3 checkm_summary.py <input.tsv> <output.txt>")

    stats = load_stats(sys.argv[1])
    write_summary(stats, sys.argv[2])


if __name__ == "__main__":
    main()

