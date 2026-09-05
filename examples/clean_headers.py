#!/usr/bin/env python3
"""Use taxutils' shared backend to clean FASTA headers."""

import argparse

from taxutils.clean_headers import clean_fasta_headers


def main():
    parser = argparse.ArgumentParser(
        description="Replace FASTA headers with cleaned accession-only headers."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA; omit it for an atomic in-place rewrite.",
    )
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()
    clean_fasta_headers(args.input, args.output, verbose=args.verbose)


if __name__ == "__main__":
    main()
