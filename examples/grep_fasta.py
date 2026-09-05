#!/usr/bin/env python3
"""Use taxutils' shared backend to select FASTA records by accession."""

import argparse

from taxutils.grep_fasta import grep_fasta


def main():
    parser = argparse.ArgumentParser(
        description="Extract FASTA records whose headers match requested accessions."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument(
        "-a",
        "--accessions",
        required=True,
        help="One accession, comma-separated accessions, or an accession file.",
    )
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--no-version", action="store_false", dest="version")
    parser.add_argument("--batch-size", type=int, default=1_000_000)
    parser.add_argument("--verbose", action="store_true")
    parser.set_defaults(version=True)
    args = parser.parse_args()

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
        f"requested={totals['requested']} scanned={totals['scanned']} "
        f"matched={totals['matched']} "
        f"missing_accession={totals['missing_accession']}"
    )


if __name__ == "__main__":
    main()
