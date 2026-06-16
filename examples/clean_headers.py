import argparse
import os
import tempfile

from taxutils import taxutils

BUFFER_SIZE = 1_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Replace FASTA headers with cleaned accession-only headers."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Path to input FASTA file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Path to output FASTA file. Defaults to overwriting the input file.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print FASTA headers where no accession is found.",
    )
    return parser.parse_args()


def first_accessions(headers, line_numbers, tu, verbose=False):
    accessions = tu.parse_accession(headers)
    for accession, header, line_number in zip(accessions, headers, line_numbers):
        if accession == "NA":
            if verbose:
                print(f"NA accession line {line_number}: {header.strip()}")
            raise ValueError(f"No accession found in FASTA header on line {line_number}: {header.strip()}")
    return accessions


def write_buffer(out_f, lines, header_positions, header_lines, header_line_numbers, tu, verbose=False):
    accessions = first_accessions(header_lines, header_line_numbers, tu, verbose=verbose)
    for position, accession in zip(header_positions, accessions):
        lines[position] = f">{accession}\n"
    out_f.writelines(lines)


def write_clean_fasta(input_path, output_path, tu, verbose=False):
    lines = []
    header_positions = []
    header_lines = []
    header_line_numbers = []
    buffer_size = 0

    def flush(out_f):
        nonlocal lines, header_positions, header_lines, header_line_numbers, buffer_size
        if lines:
            write_buffer(
                out_f,
                lines,
                header_positions,
                header_lines,
                header_line_numbers,
                tu,
                verbose=verbose,
            )
            lines = []
            header_positions = []
            header_lines = []
            header_line_numbers = []
            buffer_size = 0

    with open(input_path) as in_f, open(output_path, "w") as out_f:
        for line_number, line in enumerate(in_f, start=1):
            if line.startswith(">"):
                header_positions.append(len(lines))
                header_lines.append(line)
                header_line_numbers.append(line_number)
            lines.append(line)
            buffer_size += len(line)
            if buffer_size >= BUFFER_SIZE:
                flush(out_f)
        flush(out_f)


def overwrite_clean_fasta(input_path, tu, verbose=False):
    input_dir = os.path.dirname(os.path.abspath(input_path))
    fd, tmp_path = tempfile.mkstemp(prefix=".clean.", suffix=".fasta", dir=input_dir)
    os.close(fd)
    try:
        write_clean_fasta(input_path, tmp_path, tu, verbose=verbose)
        os.replace(tmp_path, input_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return input_path


def clean_fasta_headers(input_path, output_path=None, verbose=False):
    tu = taxutils()
    if output_path is None or os.path.abspath(output_path) == os.path.abspath(input_path):
        return overwrite_clean_fasta(input_path, tu, verbose=verbose)

    write_clean_fasta(input_path, output_path, tu, verbose=verbose)
    return output_path


def main():
    args = parse_args()
    clean_fasta_headers(args.input, args.output, verbose=args.verbose)


if __name__ == "__main__":
    main()
