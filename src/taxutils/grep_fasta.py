import argparse
import os

from taxutils.backend import call_rust


def parse_args(argv=None):
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
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Worker threads to use. Defaults to all logical CPUs.",
    )
    parser.set_defaults(version=True)
    return parser.parse_args(argv)


def grep_fasta(
    input_path,
    accession_query,
    output_path,
    version=True,
    batch_size=1_000_000,
    verbose=False,
    threads=None,
):
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")
    requested, scanned, matched, missing_accession = call_rust(
        "grep_fasta",
        os.fspath(input_path),
        os.fspath(accession_query),
        os.fspath(output_path),
        version,
        batch_size,
        verbose,
        threads,
    )
    return {
        "requested": requested,
        "scanned": scanned,
        "matched": matched,
        "missing_accession": missing_accession,
    }


def main(argv=None):
    args = parse_args(argv)
    totals = grep_fasta(
        args.input,
        args.accessions,
        args.output,
        version=args.version,
        batch_size=args.batch_size,
        verbose=args.verbose,
        threads=args.threads,
    )
    print(
        "Finished grepping FASTA: "
        f"requested={totals['requested']} scanned={totals['scanned']} "
        f"matched={totals['matched']} "
        f"missing_accession={totals['missing_accession']}"
    )


if __name__ == "__main__":
    main()
