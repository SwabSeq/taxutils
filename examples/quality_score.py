import argparse
import os

import pandas as pd

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate taxonomy distance quality metrics from Kraken output."
    )
    parser.add_argument(
        "-i",
        "--kraken-output",
        required=True,
        help="Path to Kraken standard read-level classification output.",
    )
    parser.add_argument(
        "-o",
        "--results",
        required=True,
        help="Path to output CSV.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print Kraken accession strings where no accession is found.",
    )
    return parser.parse_args()


def aggregate_lca_mapping(lca_mapping):
    counts = {}
    for item in str(lca_mapping).split():
        if ":" not in item:
            continue
        taxon, count = item.split(":", 1)
        taxon = -1 if taxon == "A" else int(taxon)
        counts[taxon] = counts.get(taxon, 0) + int(count)
    return counts


def distance_metrics(tu, labeled_taxon, kmer_counts, topology_scale):
    kmer_counts.pop(0, None)
    kmer_counts.pop(-1, None)
    classified_kmer_count = sum(kmer_counts.values())

    total_distance_moved = 0
    total_upward_loss = 0
    total_lateral_shift = 0
    distance_sq_sum = 0

    for taxon, count in kmer_counts.items():
        lca_taxon = tu.get_lca(labeled_taxon, taxon)
        upward_loss = tu.get_distance(labeled_taxon, lca_taxon)
        lateral_shift = tu.get_distance(taxon, lca_taxon)
        distance = upward_loss + lateral_shift

        total_distance_moved += count * distance
        total_upward_loss += count * upward_loss
        total_lateral_shift += count * lateral_shift
        distance_sq_sum += count * distance * distance

    average_distance_moved = (
        total_distance_moved / classified_kmer_count
        if classified_kmer_count
        else 0
    )
    average_upward_loss = (
        total_upward_loss / classified_kmer_count if classified_kmer_count else 0
    )
    average_lateral_shift = (
        total_lateral_shift / classified_kmer_count if classified_kmer_count else 0
    )
    total_normalized_distance = total_distance_moved / topology_scale
    average_normalized_distance = (
        total_normalized_distance / classified_kmer_count
        if classified_kmer_count
        else 0
    )
    distance_variance = (
        distance_sq_sum / classified_kmer_count - average_distance_moved ** 2
        if classified_kmer_count
        else 0
    )

    return {
        "total_distance_moved": total_distance_moved,
        "average_distance_moved": average_distance_moved,
        "topology_scale": topology_scale,
        "total_normalized_distance": total_normalized_distance,
        "average_normalized_distance": average_normalized_distance,
        "total_upward_loss": total_upward_loss,
        "average_upward_loss": average_upward_loss,
        "total_lateral_shift": total_lateral_shift,
        "average_lateral_shift": average_lateral_shift,
        "distance_variance": max(distance_variance, 0),
    }


def main():
    args = parse_args()
    tu = taxutils(low_memory=False)

    out = pd.read_csv(
        args.kraken_output,
        sep="\t",
        header=None,
        names=["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"],
    )

    parsed_accessions = tu.parse_accession(out["accession"])
    if args.verbose:
        for value in out.loc[parsed_accessions == "NA", "accession"]:
            print(f"NA accession: {value}")

    out["accession"] = parsed_accessions
    out = out[out["accession"] != "NA"].copy()

    tu.load_a2t(out["accession"].drop_duplicates().tolist())
    out["labeled_taxon"] = out["accession"].map(tu.a2t)
    out = out.dropna(subset=["labeled_taxon"]).astype({"labeled_taxon": int})

    topology_scales = tu.topology(
        out["labeled_taxon"].drop_duplicates(),
        anchor_rank="F",
        stat="topology_scale",
    ).to_dict()

    rows = []
    for row in out[["accession", "labeled_taxon", "lca_mapping"]].to_dict("records"):
        metrics = distance_metrics(
            tu=tu,
            labeled_taxon=row["labeled_taxon"],
            kmer_counts=aggregate_lca_mapping(row["lca_mapping"]),
            topology_scale=topology_scales[row["labeled_taxon"]],
        )
        rows.append({"accession": row["accession"], **metrics})

    results = pd.DataFrame.from_records(rows)[[
        "accession",
        "total_distance_moved",
        "average_distance_moved",
        "topology_scale",
        "total_normalized_distance",
        "average_normalized_distance",
        "total_upward_loss",
        "average_upward_loss",
        "total_lateral_shift",
        "average_lateral_shift",
        "distance_variance",
    ]]

    results_dir = os.path.dirname(args.results)
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)
    results.to_csv(args.results, index=False)

if __name__ == "__main__":
    main()
