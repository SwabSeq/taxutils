"""Remove later FASTA records sharing a parsed, versioned accession."""

import argparse
import os

from taxutils.backend import call_rust


def deduplicate_fasta(input_path, output_path=None, threads=None):
    """Keep first occurrences, returning counts of kept and removed records.

    Preserve record bytes and order. Missing accessions raise ValueError.
    Replace the destination atomically, defaulting to the input file.
    """
    kept, removed = call_rust(
        "deduplicate_fasta",
        os.fspath(input_path),
        None if output_path is None else os.fspath(output_path),
        threads,
    )
    return {"kept": kept, "removed": removed}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output FASTA. Defaults to an atomic in-place rewrite.",
    )
    parser.add_argument("--threads", type=int, default=None)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    stats = deduplicate_fasta(args.input, args.output, args.threads)
    print(f"Finished deduplicating FASTA: kept={stats['kept']} removed={stats['removed']}")


if __name__ == "__main__":
    main()
