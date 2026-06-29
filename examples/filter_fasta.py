import argparse
import os

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter records from a FASTA using accession-to-taxid lookup."
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
    taxa_group = parser.add_mutually_exclusive_group(required=True)
    taxa_group.add_argument(
        "--keep-taxids",
        help=(
            "Taxid or taxa to keep. Pass a single taxid, comma-separated taxids, "
            "or a text file containing taxids."
        ),
    )
    taxa_group.add_argument(
        "--remove-taxids",
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


def parse_taxa(value, option_name):
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
            raise ValueError(f"Invalid taxid in {option_name}: {token}") from exc

    if not taxa:
        raise ValueError(f"{option_name} did not contain any taxids")

    return taxa


def parse_taxa_to_remove(value):
    return parse_taxa(value, "--remove-taxids")


def parse_taxa_to_keep(value):
    return parse_taxa(value, "--keep-taxids")


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


def flush_records(out_f, records, tu, filter_taxa, filter_mode="remove", verbose=False):
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
            if filter_mode == "keep":
                removed += 1
                continue
            lines.append(header)
            lines.extend(sequence)
            kept += 1
            continue

        taxid = tu.a2t.get(accession)
        if taxid is None:
            missing_taxid += 1
            if verbose:
                print(f"No taxid found for {accession}: {header.strip()}")
            if filter_mode == "keep":
                removed += 1
                continue
            lines.append(header)
            lines.extend(sequence)
            kept += 1
            continue

        taxid = int(taxid)
        if filter_mode == "remove" and taxid in filter_taxa:
            removed += 1
            continue
        if filter_mode == "keep" and taxid not in filter_taxa:
            removed += 1
            continue

        lines.append(header)
        lines.extend(sequence)
        kept += 1

    out_f.writelines(lines)
    return kept, removed, missing_accession, missing_taxid


def filter_fasta(
    input_path,
    output_path,
    filter_taxa,
    filter_mode="remove",
    batch_size=5000,
    verbose=False,
):
    if filter_mode not in {"keep", "remove"}:
        raise ValueError("filter_mode must be 'keep' or 'remove'")

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
                    filter_taxa,
                    filter_mode=filter_mode,
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
                filter_taxa,
                filter_mode=filter_mode,
                verbose=verbose,
            )
            totals["kept"] += kept
            totals["removed"] += removed
            totals["missing_accession"] += missing_accession
            totals["missing_taxid"] += missing_taxid

    return totals


def main():
    args = parse_args()
    if args.keep_taxids:
        filter_taxa = parse_taxa_to_keep(args.keep_taxids)
        filter_mode = "keep"
    else:
        filter_taxa = parse_taxa_to_remove(args.remove_taxids)
        filter_mode = "remove"

    totals = filter_fasta(
        args.input,
        args.output,
        filter_taxa=filter_taxa,
        filter_mode=filter_mode,
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
