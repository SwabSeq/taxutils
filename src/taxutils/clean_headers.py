import argparse
import os

from taxutils.backend import call_rust


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Replace FASTA headers with cleaned accession-only headers."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output FASTA. Defaults to an atomic in-place rewrite.",
    )
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Worker threads to use. Defaults to all logical CPUs.",
    )
    return parser.parse_args(argv)


def clean_fasta_headers(input_path, output_path=None, verbose=False, threads=None):
    rust_output = None if output_path is None else os.fspath(output_path)
    call_rust(
        "clean_fasta_headers",
        os.fspath(input_path),
        rust_output,
        verbose,
        threads,
    )
    if output_path is None or os.path.abspath(output_path) == os.path.abspath(input_path):
        return input_path
    return output_path


def main(argv=None):
    args = parse_args(argv)
    clean_fasta_headers(args.input, args.output, verbose=args.verbose, threads=args.threads)


if __name__ == "__main__":
    main()
