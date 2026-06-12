import argparse
import os
import tempfile

from taxutils import taxutils

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
    return parser.parse_args()


def clean_header(line, tu, line_number):
    accessions = tu.parse_accession(line)
    if not accessions:
        raise ValueError(f"No accession found in FASTA header on line {line_number}: {line.strip()}")
    return f">{accessions[0]}\n"


def write_clean_fasta(input_path, output_path, tu):
    buffer = []
    buffer_size = 0

    def flush(out_f):
        nonlocal buffer, buffer_size
        if buffer:
            out_f.writelines(buffer)
            buffer = []
            buffer_size = 0

    with open(input_path) as in_f, open(output_path, "w") as out_f:
        for line_number, line in enumerate(in_f, start=1):
            if line.startswith(">"):
                line = clean_header(line, tu, line_number)
            buffer.append(line)
            buffer_size += len(line)
            if buffer_size >= 1_000_000:
                flush(out_f)
        flush(out_f)


def overwrite_clean_fasta(input_path, tu):
    input_dir = os.path.dirname(os.path.abspath(input_path))
    fd, tmp_path = tempfile.mkstemp(prefix=".clean.", suffix=".fasta", dir=input_dir)
    os.close(fd)
    try:
        write_clean_fasta(input_path, tmp_path, tu)
        os.replace(tmp_path, input_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return input_path


def clean_fasta_headers(input_path, output_path=None):
    tu = taxutils()
    if output_path is None or os.path.abspath(output_path) == os.path.abspath(input_path):
        return overwrite_clean_fasta(input_path, tu)

    write_clean_fasta(input_path, output_path, tu)
    return output_path


def main():
    args = parse_args()
    clean_fasta_headers(args.input, args.output)


if __name__ == "__main__":
    main()
