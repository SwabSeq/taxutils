"""Constants and small shared helpers.

`TAXUTILS_GLOBALS` is the one environment variable this package reads: it names
the directory holding `names.dmp`, `nodes.dmp`, `targets.json` and the accession
database. Everything else is passed as an argument.
"""

import logging
import os
import re

# Taxid reserved for reads with no assignment.
UNCLASSIFIED = 0

# Where taxonomy resources are downloaded to when a caller does not say.
DEFAULT_SAVE_FOLDER = "./taxutils/"

PATHOGEN_DICT_URLS = [
    "https://web.cs.ucla.edu/~wob/projects/taxutils/targets.json",
]

RANK_ORDER = {"U": 0, "R": 1, "D": 2, "K": 3, "P": 4, "C": 5, "O": 6, "F": 7, "G": 8, "S": 9}
CANONICAL_RANK_NAMES = {
    "U": "unclassified",
    "R": "root",
    "D": "domain",
    "K": "kingdom",
    "P": "phylum",
    "C": "class",
    "O": "order",
    "F": "family",
    "G": "genus",
    "S": "species",
}
RANK_ALIASES = {
    "U": "U",
    "UNCLASSIFIED": "U",
    "R": "R",
    "ROOT": "R",
    "D": "D",
    "DOMAIN": "D",
    "SUPERKINGDOM": "D",
    "REALM": "D",
    "K": "K",
    "KINGDOM": "K",
    "P": "P",
    "PHYLUM": "P",
    "C": "C",
    "CLASS": "C",
    "CLADE": "C",
    "O": "O",
    "ORDER": "O",
    "F": "F",
    "FAMILY": "F",
    "SUBFAMILY": "F",
    "G": "G",
    "GENUS": "G",
    "S": "S",
    "SPECIES": "S",
}
MAJOR_RANK_TO_CODE = {
    "root": "R",
    "acellular root": "R",
    "cellular root": "R",
    "domain": "D",
    "superkingdom": "D",
    "realm": "D",
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
}
ACCESSION_PATTERN = re.compile(
    r"(?<![A-Za-z0-9_])("
    r"(?:[A-Z]{2}_[A-Z]{2}[0-9]{5,}|"
    r"[A-Z]{2}_[A-Z]{1,6}[0-9]{5,}(?:[A-Z]{0,2})?|"
    r"[A-Z]{1,4}_?[0-9]{5,}|"
    r"[A-Z]{4,6}[0-9]{8,}(?:[A-Z]{0,2})?|"
    r"[A-Z]{3}[0-9]{5}|"
    r"[A-Z][0-9][A-Z0-9]{8}|"
    r"[A-Z][0-9][A-Z0-9]{3}[0-9])"
    r")(?:\.([0-9]+))?(?![A-Za-z0-9_])"
)


def default_save_folder():
    """Return the save directory named by TAXUTILS_GLOBALS, or the default.

    Read on each call rather than cached at import, so setting the variable
    after importing the package still takes effect.
    """
    return os.path.expanduser(os.environ.get("TAXUTILS_GLOBALS", DEFAULT_SAVE_FOLDER))


def resolve_save_folder(save_folder=None):
    """Return an explicit save folder, falling back to the environment."""
    if save_folder is None:
        return default_save_folder()
    return os.path.expanduser(os.fspath(save_folder))


def get_logger(name):
    logger = logging.getLogger(name)
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s : %(message)s",
        level=logging.INFO,
    )
    return logger
