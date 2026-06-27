import argparse
import os

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Remove records from a FASTA using accession-to-taxid lookup."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="/shares/swabseq/db/databases/metagenomic_data/NT_Viral_Data/AllNucleotide.fa",
        help="Path to input FASTA.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to write filtered FASTA.",
    )
    parser.add_argument(
        "--remove-taxids",
        required=True,
        help=(
            "Taxid or taxa to remove. Pass a single taxid, comma-separated taxids, "
            "or a text file containing taxids."
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5000,
        help="Number of FASTA records to resolve per lookup batch.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print headers where no accession or taxid is found.",
    )
    return parser.parse_args()


def parse_taxa_to_remove(value):
    if os.path.exists(value):
        with open(value) as f:
            text = f.read()
    else:
        text = value

    taxa = set()
    for token in text.replace(",", " ").split():
        try:
            taxa.add(int(token))
        except ValueError as exc:
            raise ValueError(f"Invalid taxid in --remove-taxids: {token}") from exc

    if not taxa:
        raise ValueError("--remove-taxids did not contain any taxids")

    return taxa


def iter_fasta_records(fasta_path):
    header = None
    sequence = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, sequence
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            yield header, sequence


def flush_records(out_f, records, tu, remove_taxa, verbose=False):
    headers = [header for header, _ in records]
    accessions = tu.parse_accession(headers)
    lookup_accessions = sorted({
        accession for accession in accessions if accession != "NA"
    })
    if lookup_accessions:
        tu.load_a2t(lookup_accessions)

    kept = 0
    removed = 0
    missing_accession = 0
    missing_taxid = 0
    lines = []

    for (header, sequence), accession in zip(records, accessions):
        if accession == "NA":
            missing_accession += 1
            if verbose:
                print(f"No accession found: {header.strip()}")
            lines.append(header)
            lines.extend(sequence)
            kept += 1
            continue

        taxid = tu.a2t.get(accession)
        if taxid is None:
            missing_taxid += 1
            if verbose:
                print(f"No taxid found for {accession}: {header.strip()}")
            lines.append(header)
            lines.extend(sequence)
            kept += 1
            continue

        if int(taxid) in remove_taxa:
            removed += 1
            continue

        lines.append(header)
        lines.extend(sequence)
        kept += 1

    out_f.writelines(lines)
    return kept, removed, missing_accession, missing_taxid


def filter_fasta(input_path, output_path, remove_taxa, batch_size=5000, verbose=False):
    tu = taxutils(low_memory=False)
    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    totals = {
        "kept": 0,
        "removed": 0,
        "missing_accession": 0,
        "missing_taxid": 0,
    }
    records = []

    with open(output_path, "w") as out_f:
        for record in iter_fasta_records(input_path):
            records.append(record)
            if len(records) >= batch_size:
                kept, removed, missing_accession, missing_taxid = flush_records(
                    out_f,
                    records,
                    tu,
                    remove_taxa,
                    verbose=verbose,
                )
                totals["kept"] += kept
                totals["removed"] += removed
                totals["missing_accession"] += missing_accession
                totals["missing_taxid"] += missing_taxid
                records = []

        if records:
            kept, removed, missing_accession, missing_taxid = flush_records(
                out_f,
                records,
                tu,
                remove_taxa,
                verbose=verbose,
            )
            totals["kept"] += kept
            totals["removed"] += removed
            totals["missing_accession"] += missing_accession
            totals["missing_taxid"] += missing_taxid

    return totals


def main():
    args = parse_args()
    remove_taxa = parse_taxa_to_remove(args.remove_taxids)
    totals = filter_fasta(
        args.input,
        args.output,
        remove_taxa=remove_taxa,
        batch_size=args.batch_size,
        verbose=args.verbose,
    )
    print(
        "Finished filtering FASTA: "
        f"kept={totals['kept']} "
        f"removed={totals['removed']} "
        f"missing_accession={totals['missing_accession']} "
        f"missing_taxid={totals['missing_taxid']}"
    )


if __name__ == "__main__":
    main()
