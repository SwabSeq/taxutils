#!/usr/bin/env python3
"""Use taxutils' shared backend to extract FASTA accessions."""

import argparse

from taxutils.extract_headers import DEFAULT_BATCH_SIZE, extract_accessions


def main():
    parser = argparse.ArgumentParser(
        description="Extract one accession per FASTA header into a single-column file."
    )
    parser.add_argument("fasta", help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    args = parser.parse_args()

    count = extract_accessions(args.fasta, args.output, batch_size=args.batch_size)
    print(f"Wrote {count} accessions to {args.output}")


if __name__ == "__main__":
    main()
