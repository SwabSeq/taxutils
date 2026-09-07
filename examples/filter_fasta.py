#!/usr/bin/env python3
"""Use taxutils' shared backend to filter FASTA records by taxid."""

import argparse

from taxutils.filter_fasta import filter_fasta, parse_taxa_to_keep, parse_taxa_to_remove


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTA records using accession-to-taxid lookup."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    taxa = parser.add_mutually_exclusive_group(required=True)
    taxa.add_argument("--keep-taxids")
    taxa.add_argument("--remove-taxids")
    parser.add_argument("--batch-size", type=int, default=5_000)
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    if args.keep_taxids:
        filter_taxa = parse_taxa_to_keep(args.keep_taxids)
        filter_mode = "keep"
    else:
        filter_taxa = parse_taxa_to_remove(args.remove_taxids)
        filter_mode = "remove"

    totals = filter_fasta(
        args.input,
        args.output,
        filter_taxa,
        filter_mode=filter_mode,
        batch_size=args.batch_size,
        verbose=args.verbose,
    )
    print(
        "Finished filtering FASTA: "
        f"kept={totals['kept']} removed={totals['removed']} "
        f"missing_accession={totals['missing_accession']} "
        f"missing_taxid={totals['missing_taxid']}"
    )


if __name__ == "__main__":
    main()
