import argparse
import os
import re

import pandas as pd

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Rewrite Kraken FASTA headers and build a Kraken2 seqid2taxid.map."
    )
    parser.add_argument(
        "-f",
        "--fasta",
        required=True,
        help="Path to input FASTA.",
    )
    parser.add_argument(
        "-k",
        "--kraken-output",
        required=True,
        help="Path to label table with accession and labeled_taxid columns.",
    )
    parser.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        help="Path to write FASTA with updated Kraken taxids.",
    )
    parser.add_argument(
        "-m",
        "--map-output",
        required=True,
        help="Path to write seqid2taxid.map.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print headers or table rows that cannot be mapped.",
    )
    return parser.parse_args()


def read_label_table(path):
    with open(path) as f:
        first_line = f.readline().strip()

    sep = "," if first_line.count(",") > first_line.count("\t") else "\t"
    fields = first_line.split(sep)
    known_columns = {
        "accession",
        "labeled_taxid",
        "labeled_taxon",
        "predicted_taxid",
        "taxon",
    }

    if known_columns.intersection(fields):
        return pd.read_csv(path, sep=sep)

    if sep == "\t" and len(fields) >= 5:
        return pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"],
        )

    raise ValueError(
        "Could not find a header with accession and labeled_taxid columns, "
        "or parse the file as standard five-column Kraken output."
    )


def find_column(columns, choices, label):
    for choice in choices:
        if choice in columns:
            return choice
    valid = ", ".join(choices)
    raise ValueError(f"Could not find {label} column. Expected one of: {valid}")


def build_accession_to_taxon(label_path, tu, verbose=False):
    labels = read_label_table(label_path)
    accession_col = find_column(labels.columns, ["accession", "seqid", "sequence_id"], "accession")
    taxon_col = find_column(
        labels.columns,
        ["labeled_taxid", "labeled_taxon", "predicted_taxid", "taxon"],
        "taxon",
    )

    labels = labels[[accession_col, taxon_col] + [
        col for col in ["total_normalized_distance"] if col in labels.columns
    ]].copy()
    labels["accession"] = tu.parse_accession(labels[accession_col])
    labels["taxon"] = pd.to_numeric(labels[taxon_col], errors="coerce")
    labels = labels[labels["accession"].ne("NA") & labels["taxon"].notna()].copy()

    if "total_normalized_distance" in labels.columns:
        labels = labels.sort_values("total_normalized_distance")
    elif verbose and labels["accession"].duplicated().any():
        print("Duplicate accessions found without total_normalized_distance; keeping first label.")

    labels = labels.drop_duplicates("accession", keep="first")
    labels["taxon"] = labels["taxon"].astype(int)
    return dict(zip(labels["accession"], labels["taxon"]))


def replace_kraken_taxid(header, taxon):
    return re.sub(
        r"(kraken:taxid\|)[0-9]+(\|)",
        rf"\g<1>{int(taxon)}\2",
        header,
        count=1,
        flags=re.IGNORECASE,
    )


def write_fixed_fasta(input_fasta, output_fasta, a2t, tu, verbose=False):
    output_dir = os.path.dirname(output_fasta)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(input_fasta) as in_f, open(output_fasta, "w") as out_f:
        for line in in_f:
            if line.startswith(">") and "kraken" in line.lower():
                accession = tu.parse_accession(line)
                taxon = a2t.get(accession)
                if taxon is None:
                    if verbose:
                        print(f"No taxid for header accession: {line.strip()}")
                else:
                    line = replace_kraken_taxid(line, taxon)
            out_f.write(line)


def write_seqid2taxid_map(fasta_path, map_output, a2t, tu, verbose=False):
    output_dir = os.path.dirname(map_output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(fasta_path) as in_f, open(map_output, "w") as out_f:
        for line in in_f:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            accession = tu.parse_accession(header)
            taxon = a2t.get(accession)
            if taxon is None:
                if verbose:
                    print(f"No taxid for map header: {header}")
                continue
            out_f.write(f"{header}\t{int(taxon)}\n")


def main():
    args = parse_args()
    tu = taxutils()
    a2t = build_accession_to_taxon(args.kraken_output, tu, verbose=args.verbose)
    write_fixed_fasta(args.fasta, args.output_fasta, a2t, tu, verbose=args.verbose)
    write_seqid2taxid_map(args.output_fasta, args.map_output, a2t, tu, verbose=args.verbose)


if __name__ == "__main__":
    main()
