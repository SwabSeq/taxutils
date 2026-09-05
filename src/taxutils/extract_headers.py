#!/usr/bin/env python3

import argparse
import os

from taxutils.backend import call_rust


DEFAULT_BATCH_SIZE = 10_000


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Extract one accession per FASTA header into a single-column file."
    )
    parser.add_argument("fasta", help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    return parser.parse_args(argv)


def extract_accessions(fasta_path, output_path, batch_size=DEFAULT_BATCH_SIZE):
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")
    return call_rust(
        "extract_accessions",
        os.fspath(fasta_path),
        os.fspath(output_path),
        batch_size,
    )


def main(argv=None):
    args = parse_args(argv)
    count = extract_accessions(args.fasta, args.output, batch_size=args.batch_size)
    print(f"Wrote {count} accessions to {args.output}")


if __name__ == "__main__":
    main()
