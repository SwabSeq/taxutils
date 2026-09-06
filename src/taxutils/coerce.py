"""Input coercion that preserves the caller's container type.

Public methods accept a scalar, list, ``pandas.Series`` or ``numpy.ndarray``
and return the same shape. The per-method conversions differ in dtype handling
(``astype(bool)`` versus ``np.asarray(...).reshape``), so only the genuinely
shared helpers live here; each method keeps its own return conversion.
"""

import numpy as np
import pandas as pd


def as_taxa_list(taxa):
    """Return taxa as a list of ints, dropping missing values."""
    if np.isscalar(taxa):
        values = [taxa]
    elif isinstance(taxa, pd.Series):
        values = taxa.tolist()
    elif isinstance(taxa, np.ndarray):
        values = taxa.ravel().tolist()
    else:
        values = list(taxa)
    return [int(t) for t in values if not pd.isna(t)]


def is_taxon_scalar(value):
    return value is None or value is pd.NA or np.isscalar(value)


def pairwise_values(value):
    if isinstance(value, pd.Series):
        return value.tolist()
    if isinstance(value, np.ndarray):
        return value.ravel().tolist()
    return list(value)


def pairwise_taxa_result(taxon_a, taxon_b, check, name):
    """Apply `check` across two taxon inputs, preserving the input container."""
    scalar_a = is_taxon_scalar(taxon_a)
    scalar_b = is_taxon_scalar(taxon_b)

    if scalar_a and scalar_b:
        return check(taxon_a, taxon_b)

    if scalar_a != scalar_b:
        raise ValueError("taxon_a and taxon_b must both be scalar or both be list-like")

    values_a = pairwise_values(taxon_a)
    values_b = pairwise_values(taxon_b)
    if len(values_a) != len(values_b):
        raise ValueError("taxon_a and taxon_b must have the same length")

    values = [check(a, b) for a, b in zip(values_a, values_b)]

    if isinstance(taxon_a, pd.Series):
        return pd.Series(values, index=taxon_a.index, name=name, dtype=bool)

    if isinstance(taxon_a, np.ndarray):
        return np.asarray(values, dtype=bool).reshape(taxon_a.shape)

    return values
