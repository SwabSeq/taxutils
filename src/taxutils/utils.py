# global_utils.py

import os
import re

TAXUTILS_GLOBALS = dict()
TAXUTILS_GLOBALS["save_folder"] = os.path.expanduser(os.environ.get("TAXUTILS_GLOBALS", "./taxutils/"))
TAXUTILS_GLOBALS["pathogen_dict_urls"] = [
    "https://web.cs.ucla.edu/~wob/projects/taxutils/targets.json",
]
TAXUTILS_GLOBALS["UNCLASSIFIED"] = 0

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

def get_logger(name):
    import logging
    logger = logging.getLogger(name)
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s : %(message)s",
        level=logging.INFO,
    )
    return logger
