import argparse
import os

import pandas as pd

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Enumerate candidate label_taxon values from Kraken LCA mappings."
    )
    parser.add_argument(
        "-i",
        "--kraken-output",
        required=True,
        help="Path to Kraken standard read-level classification output.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to write enumerated label table.",
    )
    parser.add_argument(
        "--exclude-original",
        action="store_true",
        help="Only write enumerated candidate-label rows.",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum kmer count required for a taxon to be enumerated.",
    )
    return parser.parse_args()


def candidate_taxa(lca_mapping, min_count=1):
    taxa = []
    seen = set()
    for item in str(lca_mapping).split():
        if ":" not in item:
            continue
        taxon, count = item.split(":", 1)
        if taxon in {"0", "A"}:
            continue
        taxon = int(taxon)
        if taxon in seen or int(count) < min_count:
            continue
        seen.add(taxon)
        taxa.append(taxon)
    return taxa


def enumerate_labels(df, keep_original=True, min_count=1):
    tu = taxutils(low_memory=False, keep_accession_downloads=False)
    df = df.copy()
    df["accession"] = tu.parse_accession(df["accession"])
    df = df[df["accession"] != "NA"].copy()
    tu.load_a2t(df["accession"].drop_duplicates().tolist())
    df["a2t_taxon"] = pd.to_numeric(df["accession"].map(tu.a2t), errors="coerce")
    df = df.dropna(subset=["a2t_taxon"]).astype({"a2t_taxon": int})

    rows = []
    for row in df.to_dict("records"):
        wrote_lca_mapping = False
        if keep_original:
            original = row.copy()
            original["label_taxon"] = original["a2t_taxon"]
            rows.append(original)
            wrote_lca_mapping = True
        for taxon in candidate_taxa(row["lca_mapping"], min_count=min_count):
            enumerated = row.copy()
            enumerated["label_taxon"] = taxon
            if wrote_lca_mapping:
                enumerated["lca_mapping"] = pd.NA
            else:
                wrote_lca_mapping = True
            rows.append(enumerated)
    return pd.DataFrame.from_records(rows)[
        ["status", "accession", "label_taxon", "seqlen", "lca_mapping", "a2t_taxon"]
    ]


def main():
    args = parse_args()
    df = pd.read_csv(
        args.kraken_output,
        sep="\t",
        header=None,
        names=["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"],
    )
    enumerated = enumerate_labels(
        df,
        keep_original=not args.exclude_original,
        min_count=args.min_count,
    )

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    enumerated.to_csv(args.output, sep="\t", header=False, index=False)


if __name__ == "__main__":
    main()
