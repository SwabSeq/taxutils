"""Downloading and loading taxonomy resources, and accession lookups.

Every function takes the save folder and thread count as arguments. The only
implicit value is the `TAXUTILS_GLOBALS` default resolved by
`utils.resolve_save_folder`, and that happens once at the caller's boundary.
"""

import json
import os
import shutil
import tarfile
import urllib.request

import pandas as pd

from .backend import call_rust
from .parse import accessions_for_lookup
from .ranks import assign_rank_codes
from .tree import build_children, build_parent, get_parents, get_subtree, taxonomic_order
from .utils import (
    CANONICAL_RANK_NAMES,
    PATHOGEN_DICT_URLS,
    RANK_ORDER,
    UNCLASSIFIED,
    get_logger,
    resolve_save_folder,
)

logger = get_logger(__name__)

TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"


def download_file(url, path):
    """Download `url` to `path`, leaving no partial file behind on failure."""
    tmp_path = f"{path}.tmp"
    try:
        urllib.request.urlretrieve(url, tmp_path)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def download_taxdump(save_path, names_path, nodes_path):
    """Fetch names.dmp and nodes.dmp out of the NCBI taxdump tarball."""
    logger.info(f"Downloading {names_path}, {nodes_path}...")
    tarball_path = os.path.join(save_path, "taxdump.tar.gz")
    download_file(TAXDUMP_URL, tarball_path)

    with tarfile.open(tarball_path, "r:gz") as tar:
        members = {member.name: member for member in tar.getmembers()}
        for filename, output_path in [
            ("names.dmp", names_path),
            ("nodes.dmp", nodes_path),
        ]:
            if filename not in members:
                raise RuntimeError(f"Could not find {filename} in taxdump.")
            source = tar.extractfile(members[filename])
            if source is None:
                raise RuntimeError(f"Could not extract {filename} from taxdump.")
            with source, open(output_path, "wb") as out:
                shutil.copyfileobj(source, out)

    os.remove(tarball_path)


def download_targets(targets_json):
    """Fetch the pathogen target list, trying each configured mirror."""
    for url in PATHOGEN_DICT_URLS:
        try:
            logger.info(f"Downloading targets.json from {url}...")
            download_file(url, targets_json)
            return
        except Exception as e:
            logger.warning(f"Failed to download from {url}: {e}")
    raise RuntimeError("Could not download targets.json from any URL.")


def build_names(names_path):
    names = {}
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxon = parts[0]
                name = parts[1]
                names[int(taxon)] = name
    names[UNCLASSIFIED] = "unclassified"
    return names


def build_nodes(nodes_path, names):
    nodes = pd.read_csv(
        nodes_path, sep="|", header=None, usecols=[0, 1, 2],
        names=["taxon", "parent", "rank"], dtype={"taxon": int, "parent": int, "rank": str},
        engine="python"
    )
    nodes["rank"] = nodes["rank"].str.strip().str.lower()
    parent = dict(zip(nodes["taxon"], nodes["parent"]))
    rank_map = dict(zip(nodes["taxon"], nodes["rank"]))
    rank_codes = assign_rank_codes(parent, rank_map)

    nodes["rank_code"] = nodes["taxon"].map(rank_codes)
    nodes["rank_base"] = nodes["rank_code"].str[0]
    nodes["rank_idx"] = nodes["rank_base"].map(RANK_ORDER)
    nodes["new_rank"] = nodes["rank_base"].map(CANONICAL_RANK_NAMES)

    return nodes


def build_target_taxa(nodes, names, targets_json):
    with open(targets_json) as f:
        pdict = json.load(f)
    pathogen_taxa = {int(v) for v in pdict["pathogens"].values()}

    parent = build_parent(nodes)
    tree = build_children(parent)
    rank = dict(zip(nodes["taxon"], nodes["rank_code"]))
    rank_idx = dict(zip(nodes["taxon"], nodes["rank_idx"]))

    taxa = set()
    for taxon in pathogen_taxa:
        taxa.update(get_subtree(taxon, tree))
        taxa.update(get_parents(taxon, parent, rank_idx, rank="F"))

    return taxonomic_order(taxa, parent, rank, names)


def ensure_a2t_db(
    save_folder=None,
    wgs=False,
    refresh=False,
    threads=None,
):
    """Build, refresh or validate the SQLite accession index; return its path.

    A missing or unusable database is built from scratch; `refresh` asks NCBI
    whether the sources changed and applies only the difference.
    """
    return os.fspath(
        call_rust(
            "ensure_accession_database",
            os.fspath(resolve_save_folder(save_folder)),
            wgs,
            refresh,
            threads,
        )
    )


def build_a2t(
    accessions,
    save_folder=None,
    low_memory=True,
    wgs=False,
    threads=None,
):
    """Return an accession-to-taxid map for the requested accessions."""
    return call_rust(
        "lookup_accession_taxids",
        os.fspath(resolve_save_folder(save_folder)),
        accessions_for_lookup(accessions),
        low_memory,
        wgs,
        threads,
    )


def get_t2a(
    taxa,
    save_folder=None,
    low_memory=True,
    wgs=False,
    threads=None,
):
    """Return the set of accessions assigned directly to the given taxa.

    Rust scans the compressed source directly in low-memory mode and uses the
    indexed local SQLite database otherwise.
    """
    return call_rust(
        "lookup_taxid_accessions",
        os.fspath(resolve_save_folder(save_folder)),
        [int(taxon) for taxon in taxa],
        low_memory,
        wgs,
        threads,
    )
