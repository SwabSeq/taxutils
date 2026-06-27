import argparse
import os

from taxutils import taxutils


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract FASTA records whose headers match requested accessions."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input FASTA file.",
    )
    parser.add_argument(
        "-a",
        "--accessions",
        required=True,
        help=(
            "Accession query. Pass a single accession, comma-separated accessions, "
            "or a text file containing accessions."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to output FASTA file containing matched records.",
    )
    parser.add_argument(
        "--no-version",
        action="store_false",
        dest="version",
        help="Strip accession versions before matching.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1_000_000,
        help="Approximate number of FASTA bytes to process per accession-parsing batch.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print FASTA headers where no accession is found.",
    )
    parser.set_defaults(version=True)
    return parser.parse_args()


def read_accession_query(value):
    if os.path.exists(value):
        with open(value) as f:
            text = f.read()
    else:
        text = value
    return [token for token in text.replace(",", " ").split() if token]


def parse_requested_accessions(value, tu, version=True):
    query_values = read_accession_query(value)
    accessions = {
        accession
        for accession in tu.parse_accession(query_values, version=version)
        if accession != "NA"
    }
    if not accessions:
        raise ValueError("No accessions were found in --accessions")
    return accessions


def iter_fasta_records(fasta_path, buffer_size):
    header = None
    sequence = []
    with open(fasta_path, "rb", buffering=buffer_size) as f:
        for line in f:
            if line.startswith(b">"):
                if header is not None:
                    yield header, sequence
                header = line
                sequence = []
            else:
                sequence.append(line)
        if header is not None:
            yield header, sequence


def write_matching_records(out_f, records, requested_accessions, tu, version=True, verbose=False):
    headers = [header.decode("utf-8", errors="ignore") for header, _ in records]
    accessions = tu.parse_accession(headers, version=version)
    scanned = 0
    matched = 0
    missing_accession = 0
    lines = []

    for (header, sequence), accession in zip(records, accessions):
        scanned += 1
        if accession == "NA":
            missing_accession += 1
            if verbose:
                print(f"NA accession: {header.decode('utf-8', errors='ignore').strip()}")
            continue
        if accession in requested_accessions:
            lines.append(header)
            lines.extend(sequence)
            matched += 1

    if lines:
        out_f.writelines(lines)

    return scanned, matched, missing_accession


def grep_fasta(input_path, accession_query, output_path, version=True, batch_size=1_000_000, verbose=False):
    tu = taxutils()
    requested_accessions = parse_requested_accessions(accession_query, tu, version=version)

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    totals = {
        "requested": len(requested_accessions),
        "scanned": 0,
        "matched": 0,
        "missing_accession": 0,
    }
    records = []
    buffered_bytes = 0

    with open(output_path, "wb", buffering=batch_size) as out_f:
        for header, sequence in iter_fasta_records(input_path, batch_size):
            records.append((header, sequence))
            buffered_bytes += len(header) + sum(len(line) for line in sequence)
            if buffered_bytes >= batch_size:
                scanned, matched, missing_accession = write_matching_records(
                    out_f,
                    records,
                    requested_accessions,
                    tu,
                    version=version,
                    verbose=verbose,
                )
                totals["scanned"] += scanned
                totals["matched"] += matched
                totals["missing_accession"] += missing_accession
                records = []
                buffered_bytes = 0

        if records:
            scanned, matched, missing_accession = write_matching_records(
                out_f,
                records,
                requested_accessions,
                tu,
                version=version,
                verbose=verbose,
            )
            totals["scanned"] += scanned
            totals["matched"] += matched
            totals["missing_accession"] += missing_accession

    return totals


def main():
    args = parse_args()
    totals = grep_fasta(
        args.input,
        args.accessions,
        args.output,
        version=args.version,
        batch_size=args.batch_size,
        verbose=args.verbose,
    )
    print(
        "Finished grepping FASTA: "
        f"requested={totals['requested']} "
        f"scanned={totals['scanned']} "
        f"matched={totals['matched']} "
        f"missing_accession={totals['missing_accession']}"
    )


if __name__ == "__main__":
    main()
