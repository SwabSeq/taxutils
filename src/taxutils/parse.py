"""Accession parsing from FASTA and Kraken-style headers."""

import numpy as np
import pandas as pd

from .utils import ACCESSION_PATTERN


def parse_accession(strings, version: bool = True):
    """Return the first accession parsed from each string, preserving input type.

    Strings with no recognizable accession yield ``"NA"`` rather than being
    dropped, so the result always lines up with the input.
    """
    def parse_text(text):
        if pd.isna(text):
            return "NA"
        match = ACCESSION_PATTERN.search(str(text).upper())
        if match is None:
            return "NA"
        accession = match.group(1)
        accession_version = match.group(2)
        if version and accession_version is not None:
            accession = f"{accession}.{accession_version}"
        return accession

    if strings is None or isinstance(strings, str):
        return parse_text(strings)

    if isinstance(strings, pd.Series):
        return strings.map(parse_text)

    if isinstance(strings, np.ndarray):
        values = strings.astype(object)
        return np.vectorize(parse_text, otypes=[object])(values)

    if np.isscalar(strings):
        return parse_text(strings)

    return [parse_text(text) for text in strings]


def accessions_for_lookup(accessions):
    """Return versioned accessions suitable for a backend lookup.

    Unparseable entries are dropped here because a lookup has nothing to do
    with them.
    """
    parsed = parse_accession(accessions, version=True)
    if isinstance(parsed, str):
        values = [parsed]
    elif isinstance(parsed, pd.Series):
        values = parsed.tolist()
    elif isinstance(parsed, np.ndarray):
        values = parsed.ravel().tolist()
    else:
        values = list(parsed)
    return [accession for accession in values if accession != "NA" and not pd.isna(accession)]
