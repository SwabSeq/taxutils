import argparse
import os

import pandas as pd

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze Kraken read-level classifications with taxutils taxonomy."
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
        default="./results/scan.csv",
        help="Path to results CSV.",
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


def main():
    args = parse_args()
    tu = taxutils(low_memory=False)

    print("Processing Kraken output file...")
    out = pd.read_csv(
        args.kraken_output,
        sep="\t",
        header=None,
        names=["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"],
    )
    out["accession"] = tu.parse_accession(out["accession"])
    out = out.dropna(subset=["accession"]).copy()

    tu.load_a2t(out["accession"].dropna().unique().tolist())
    out["labeled_taxid"] = out["accession"].map(tu.a2t)
    out = (
        out.dropna(subset=["labeled_taxid"])
        .astype({"labeled_taxid": int, "predicted_taxid": int})
        .assign(
            predicted_name=lambda df: df["predicted_taxid"].map(tu.names),
            labeled_name=lambda df: df["labeled_taxid"].map(tu.names),
        )
        .dropna(subset=["labeled_name"])
        .copy()
    )

    ranks = tu.nodes.set_index("taxon")["rank_code"].to_dict()
    rank_idx = tu.nodes.set_index("taxon")["rank_idx"].to_dict()
    rank_order = {rank: i for i, rank in enumerate(tu.get_rank_order())}
    out["labeled_rank"] = out["labeled_taxid"].map(ranks)

    targets = set(tu.target_taxa)
    out["is_target"] = out["labeled_taxid"].isin(targets)

    print("Scanning read-level LCA mappings...")
    metrics = []
    for row in out.to_dict("records"):
        kmer_counts = aggregate_lca_mapping(row["lca_mapping"])
        total_kmer_count = sum(kmer_counts.values())
        unclassified_kmer_count = kmer_counts.pop(0, 0)
        masked_kmers = kmer_counts.pop(-1, 0)
        classified_kmer_count = sum(kmer_counts.values())

        target_loss_count = 0
        total_distance_moved = 0
        distance_sq_sum = 0
        common_kmer_count_F = 0
        common_kmer_count_G = 0
        common_kmer_count_S = 0
        root_count = 0
        lost_details = []

        labeled_taxon = row["labeled_taxid"]
        for taxon, count in kmer_counts.items():
            lca_taxon = tu.get_lca(labeled_taxon, taxon)
            distance = tu.get_distance(labeled_taxon, taxon)
            total_distance_moved += count * distance
            distance_sq_sum += count * distance * distance

            if taxon in targets and lca_taxon not in targets:
                target_loss_count += count
                lost_details.append(
                    f"{taxon}({tu.names.get(taxon, '')}):{count}->{lca_taxon}({tu.names.get(lca_taxon, '')})"
                )

            if rank_idx.get(lca_taxon, rank_order["F"]) < rank_order["F"]:
                common_kmer_count_F += count
            if rank_idx.get(lca_taxon, rank_order["G"]) < rank_order["G"]:
                common_kmer_count_G += count
            if rank_idx.get(lca_taxon, rank_order["S"]) < rank_order["S"]:
                common_kmer_count_S += count
            if lca_taxon == 1:
                root_count += count

        average_distance_moved = (
            total_distance_moved / classified_kmer_count
            if classified_kmer_count
            else 0
        )
        distance_variance = (
            distance_sq_sum / classified_kmer_count - average_distance_moved ** 2
            if classified_kmer_count
            else 0
        )

        metrics.append({
            "total_kmer_count": total_kmer_count,
            "unclassified_kmer_count": unclassified_kmer_count,
            "masked_kmers": masked_kmers,
            "unique_taxid_count": len(kmer_counts),
            "classified_kmer_count": classified_kmer_count,
            "total_distance_moved": total_distance_moved,
            "average_distance_moved": average_distance_moved,
            "distance_variance": max(distance_variance, 0),
            "target_loss_count": target_loss_count,
            "common_kmer_count_F": common_kmer_count_F,
            "common_kmer_count_G": common_kmer_count_G,
            "common_kmer_count_S": common_kmer_count_S,
            "root_count": root_count,
            "human_kmer_count": kmer_counts.get(9606, 0),
            "lost_details": " ".join(lost_details),
        })

    out = out.join(pd.DataFrame.from_records(metrics, index=out.index))

    total = out["total_kmer_count"].replace(0, float("nan"))
    out["unclassified_ratio"] = out["unclassified_kmer_count"] / total
    out["masked_ratio"] = out["masked_kmers"] / total
    out["root_ratio"] = out["root_count"] / total
    out["F_extra_ratio"] = (out["common_kmer_count_F"] - out["root_count"]) / total
    out["G_extra_ratio"] = (out["common_kmer_count_G"] - out["common_kmer_count_F"]) / total
    out["S_extra_ratio"] = (out["common_kmer_count_S"] - out["common_kmer_count_G"]) / total
    out["remaining_ratio"] = (
        out["total_kmer_count"]
        - out["unclassified_kmer_count"]
        - out["common_kmer_count_S"]
        - out["masked_kmers"]
    ) / total
    out["target_loss_ratio"] = out["target_loss_count"] / total
    out["added_value_ratio"] = out["remaining_ratio"] + out["unclassified_ratio"]
    ratio_cols = [col for col in out.columns if col.endswith("_ratio")]
    out[ratio_cols] = out[ratio_cols].fillna(0)

    results_dir = os.path.dirname(args.results)
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)

    results_basename = os.path.splitext(os.path.basename(args.results))[0]
    accessions_path = os.path.join(results_dir, f"{results_basename}_tpcf_accessions.txt")

    print(f"Saving results to {args.results}...")
    out.to_csv(args.results, index=False)


if __name__ == "__main__":
    main()
