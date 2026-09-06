import argparse
import os

from taxutils.backend import call_rust
from taxutils.utils import resolve_save_folder


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter records from a FASTA using accession-to-taxid lookup."
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    taxa_group = parser.add_mutually_exclusive_group(required=True)
    taxa_group.add_argument("--keep-taxids")
    taxa_group.add_argument("--remove-taxids")
    parser.add_argument("--batch-size", type=int, default=5_000)
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--save-folder",
        default=None,
        help="Taxonomy resource directory. Defaults to $TAXUTILS_GLOBALS.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Worker threads to use. Defaults to all logical CPUs.",
    )
    return parser.parse_args(argv)


def parse_taxa(value, option_name):
    if os.path.exists(value):
        with open(value) as f:
            text = f.read()
    else:
        text = value

    taxa = set()
    for token in text.replace(",", " ").split():
        try:
            taxa.add(int(token))
        except ValueError as exc:
            raise ValueError(f"Invalid taxid in {option_name}: {token}") from exc
    if not taxa:
        raise ValueError(f"{option_name} did not contain any taxids")
    return taxa


def parse_taxa_to_remove(value):
    return parse_taxa(value, "--remove-taxids")


def parse_taxa_to_keep(value):
    return parse_taxa(value, "--keep-taxids")


def filter_fasta(
    input_path,
    output_path,
    filter_taxa,
    filter_mode="remove",
    batch_size=5_000,
    verbose=False,
    save_folder=None,
    threads=None,
):
    if filter_mode not in {"keep", "remove"}:
        raise ValueError("filter_mode must be 'keep' or 'remove'")
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")

    kept, removed, missing_accession, missing_taxid = call_rust(
        "filter_fasta",
        os.fspath(input_path),
        os.fspath(output_path),
        [int(taxon) for taxon in filter_taxa],
        os.fspath(resolve_save_folder(save_folder)),
        filter_mode,
        batch_size,
        verbose,
        False,
        threads,
    )
    return {
        "kept": kept,
        "removed": removed,
        "missing_accession": missing_accession,
        "missing_taxid": missing_taxid,
    }


def main(argv=None):
    args = parse_args(argv)
    if args.keep_taxids:
        filter_taxa = parse_taxa_to_keep(args.keep_taxids)
        filter_mode = "keep"
    else:
        filter_taxa = parse_taxa_to_remove(args.remove_taxids)
        filter_mode = "remove"

    totals = filter_fasta(
        args.input,
        args.output,
        filter_taxa=filter_taxa,
        filter_mode=filter_mode,
        batch_size=args.batch_size,
        verbose=args.verbose,
        save_folder=args.save_folder,
        threads=args.threads,
    )
    print(
        "Finished filtering FASTA: "
        f"kept={totals['kept']} removed={totals['removed']} "
        f"missing_accession={totals['missing_accession']} "
        f"missing_taxid={totals['missing_taxid']}"
    )


if __name__ == "__main__":
    main()
