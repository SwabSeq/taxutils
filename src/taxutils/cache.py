"""Compressed, atomic caching of built taxonomy objects."""

import gzip
import os
import pickle
import tempfile

import pandas as pd

from .utils import get_logger

logger = get_logger(__name__)

# Increment when taxonomy construction, object layout, or rank rules change.
TAXONOMY_CACHE_VERSION = 1
NODE_COLUMNS = (
    "taxon", "parent", "rank", "rank_code", "rank_base", "rank_idx", "new_rank",
)


def taxonomy_sources(*paths):
    """Fingerprint all taxonomy inputs without reading large source files."""
    sources = []
    for path in paths:
        stat = os.stat(path)
        sources.append((os.path.realpath(path), stat.st_size, stat.st_mtime_ns))
    return tuple(sources)


def load_taxonomy_cache(cache_path, sources, object_type):
    """Return a compatible cached taxonomy object, otherwise None."""
    try:
        with gzip.open(cache_path, "rb") as cached:
            header = pickle.load(cached)
            if header != (TAXONOMY_CACHE_VERSION, sources):
                return None
            tu = pickle.load(cached)
            # Consume the trailer so truncated or checksum-invalid gzip files
            # cannot be accepted even when the pickle was fully decoded.
            if cached.read(1):
                return None
        if not (
            isinstance(tu, object_type)
            and isinstance(tu.nodes, pd.DataFrame)
            and tuple(tu.nodes.columns) == NODE_COLUMNS
            and isinstance(tu.names, dict)
            and isinstance(tu.parent, dict)
            and isinstance(tu.target_taxa, list)
            and tu.a2t is None
            and tu._tree is None
            and tu._descendant_index is None
            and tu._depth == {}
            and tu._a2t_checked is False
        ):
            return None
        logger.info(f"Loading built taxonomy from {cache_path}...")
        return tu
    except FileNotFoundError:
        return None
    except Exception as error:
        logger.warning(f"Could not read taxonomy cache {cache_path}: {error}")
        return None


def save_taxonomy_cache(cache_path, sources, tu):
    """Save a built object with fast gzip compression and atomic installation."""
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=os.path.dirname(os.path.abspath(cache_path)),
            prefix=".taxutils-", suffix=".tmp", delete=False,
        ) as output:
            temporary_path = output.name
            with gzip.GzipFile(fileobj=output, mode="wb", compresslevel=1, mtime=0) as cached:
                pickle.dump((TAXONOMY_CACHE_VERSION, sources), cached, protocol=pickle.HIGHEST_PROTOCOL)
                pickle.dump(tu, cached, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary_path, cache_path)
    except Exception as error:
        logger.warning(f"Could not save taxonomy cache {cache_path}: {error}")
    finally:
        if temporary_path is not None:
            try:
                os.remove(temporary_path)
            except FileNotFoundError:
                pass
            except OSError as error:
                logger.warning(f"Could not remove temporary taxonomy cache: {error}")
