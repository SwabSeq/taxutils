# global_utils.py

import json
import os
import sys
import time
from contextlib import ContextDecorator
from datetime import datetime, timezone

try:
    import resource
except ImportError:  # pragma: no cover
    resource = None

TAXUTILS_GLOBALS = dict()
TAXUTILS["save_folder"] = os.path.expanduser(os.environ.get("TAXUTILS_GLOBALS", "./taxutils/"))
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
