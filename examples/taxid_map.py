import argparse
from pathlib import Path
import tempfile

from taxutils import taxutils
from taxutils.extract_headers import extract_accessions as extract_fasta_accessions


BUFFER_SIZE = 1_000_000


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create an accession-to-taxon map from FASTA headers."
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
        required=True,
        help="Path to output accession-to-taxon map.",
    )
    return parser.parse_args()


def extract_accessions(fasta_path):
    with tempfile.TemporaryDirectory() as temporary:
        extracted = Path(temporary) / "accessions.txt"
        extract_fasta_accessions(fasta_path, extracted)
        return list(dict.fromkeys(extracted.read_text().splitlines()))


def write_taxid_map(accessions, output_path, tu):
    tu.load_a2t(accessions)
    buffer = []
    buffer_size = 0

    def flush(out_f):
        nonlocal buffer, buffer_size
        if buffer:
            out_f.writelines(buffer)
            buffer = []
            buffer_size = 0

    with open(output_path, "w", buffering=BUFFER_SIZE) as f:
        for accession in accessions:
            taxon = tu.a2t.get(accession)
            if taxon is not None:
                line = f"{accession}\t{taxon}\n"
                buffer.append(line)
                buffer_size += len(line)
                if buffer_size >= BUFFER_SIZE:
                    flush(f)
        flush(f)


def build_taxid_map(input_path, output_path):
    tu = taxutils(low_memory=False, keep_accession_downloads=False)
    accessions = extract_accessions(input_path)
    write_taxid_map(accessions, output_path, tu)
    return output_path


def main():
    args = parse_args()
    build_taxid_map(args.input, args.output)


if __name__ == "__main__":
    main()
