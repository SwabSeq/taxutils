import argparse

from taxutils import taxutils


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


def iter_header_lines(fasta_path):
    with open(fasta_path, "rb", buffering=BUFFER_SIZE) as f:
        for line in f:
            if line.startswith(b">"):
                yield line.decode("utf-8", errors="ignore")


def extract_accessions(fasta_path, tu):
    accessions = {}
    headers = []
    headers_size = 0

    def add_headers():
        nonlocal headers, headers_size
        if headers:
            for accession in tu.parse_accession(headers):
                accessions.setdefault(accession, None)
            headers = []
            headers_size = 0

    for header in iter_header_lines(fasta_path):
        headers.append(header)
        headers_size += len(header)
        if headers_size >= BUFFER_SIZE:
            add_headers()
    add_headers()
    return list(accessions)


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
    tu = taxutils(low_memory=False)
    accessions = extract_accessions(input_path, tu)
    write_taxid_map(accessions, output_path, tu)
    return output_path


def main():
    args = parse_args()
    build_taxid_map(args.input, args.output)


if __name__ == "__main__":
    main()
