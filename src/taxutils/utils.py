# global_utils.py

import os

TAXUTILS_GLOBALS = dict()
TAXUTILS_GLOBALS["save_folder"] = os.path.expanduser(os.environ.get("TAXUTILS_GLOBALS", "./taxutils/"))
TAXUTILS_GLOBALS["pathogen_dict_urls"] = [
    "https://web.cs.ucla.edu/~wob/projects/trident/targets.json",
]
TAXUTILS_GLOBALS["UNCLASSIFIED"] = 0
TAXUTILS_GLOBALS["UNMAPPED"] = -1

def get_logger(name):
    import logging
    logger = logging.getLogger(name)
    logging.basicConfig(
        format="%(asctime)s | %(levelname)s : %(message)s",
        level=logging.INFO,
    )
    return logger
