"""Command-line interface for taxutils FASTA utilities."""

import argparse
import importlib
import sys


COMMANDS = {
    "extract": (
        "Extract one accession per FASTA header.",
        "taxutils.extract_headers",
    ),
    "clean": (
        "Replace FASTA headers with accession-only headers.",
        "taxutils.clean_headers",
    ),
    "grep": (
        "Extract FASTA records matching requested accessions.",
        "taxutils.grep_fasta",
    ),
    "filter": (
        "Filter FASTA records using accession-to-taxid lookup.",
        "taxutils.filter_fasta",
    ),
}


def build_parser():
    command_help = "\n".join(
        f"  {name:<8}{description}"
        for name, (description, _) in COMMANDS.items()
    )
    parser = argparse.ArgumentParser(
        prog="taxutils",
        description="Utilities for working with taxonomy data and FASTA files.",
        epilog=f"commands:\n{command_help}",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage="taxutils [-h] {extract,clean,grep,filter} ...",
    )
    parser.add_argument(
        "command",
        nargs="?",
        choices=COMMANDS,
        help="Command to run (omit for this help message).",
    )
    return parser


def main(argv=None):
    """Run a taxutils subcommand."""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()

    if not argv or argv[0] in {"-h", "--help"}:
        parser.print_help()
        return 0

    command = argv.pop(0)
    if command not in COMMANDS:
        parser.error(
            f"invalid choice: {command!r} "
            f"(choose from {', '.join(COMMANDS)})"
        )

    _, module_name = COMMANDS[command]
    command_main = importlib.import_module(module_name).main
    command_main(argv or ["--help"])
    return 0


if __name__ == "__main__":
    main()
