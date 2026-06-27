#!/usr/bin/env python3

import argparse
import os
import tempfile
from pathlib import Path

from taxutils.taxutils import parse_accession


DEFAULT_BATCH_SIZE = 10_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract one accession per FASTA header into a single-column file."
    )
    parser.add_argument("fasta", help="Path to the input FASTA file.")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Path to the output accession file.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f"Number of headers to parse per batch (default: {DEFAULT_BATCH_SIZE}).",
    )
    return parser.parse_args()


def write_accession_batch(out_f, headers, line_numbers):
    accessions = parse_accession(headers)
    for accession, header, line_number in zip(accessions, headers, line_numbers):
        if accession == "NA":
            raise ValueError(
                f"No accession found in FASTA header on line {line_number}: "
                f"{header.strip()}"
            )
        out_f.write(f"{accession}\n")


def extract_accessions(fasta_path, output_path, batch_size=DEFAULT_BATCH_SIZE):
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")

    fasta_path = Path(fasta_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fd, temporary_path = tempfile.mkstemp(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
        text=True,
    )
    os.close(fd)

    header_count = 0
    headers = []
    line_numbers = []

    try:
        with fasta_path.open() as in_f, open(temporary_path, "w") as out_f:
            for line_number, line in enumerate(in_f, start=1):
                if not line.startswith(">"):
                    continue

                headers.append(line)
                line_numbers.append(line_number)
                if len(headers) >= batch_size:
                    write_accession_batch(out_f, headers, line_numbers)
                    header_count += len(headers)
                    headers.clear()
                    line_numbers.clear()

            if headers:
                write_accession_batch(out_f, headers, line_numbers)
                header_count += len(headers)

        os.replace(temporary_path, output_path)
    finally:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)

    return header_count


def main():
    args = parse_args()
    count = extract_accessions(args.fasta, args.output, batch_size=args.batch_size)
    print(f"Wrote {count} accessions to {args.output}")


if __name__ == "__main__":
    main()
