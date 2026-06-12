# taxutils.py

import pandas as pd
import numpy as np
from collections import defaultdict
from typing import List
from dataclasses import dataclass
import os, json, urllib.request, tarfile, gzip
import shutil
import sqlite3

from .utils import (
    ACCESSION_PATTERN,
    CANONICAL_RANK_NAMES,
    MAJOR_RANK_TO_CODE,
    RANK_ALIASES,
    RANK_ORDER,
    TAXUTILS_GLOBALS,
    get_logger,
)

logger = get_logger(__name__)


def _download_file(url, path):
    tmp_path = f"{path}.tmp"
    try:
        urllib.request.urlretrieve(url, tmp_path)
        os.replace(tmp_path, path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

@dataclass
class TaxonomicUtils:
    names: dict
    nodes: dict
    target_taxa: list
    a2t: dict = None
    parent: dict = None
    _low_memory: bool = True

    def __post_init__(self):
        if self.parent is None:
            self.parent = build_parent(self.nodes)
        self._tree = None
        self._a2t_checked = False
    
    def __repr__(self):
        fields = []
        for f in self.__dataclass_fields__:
            if f.startswith("_"):
                continue
            val = getattr(self, f)
            if val is None:
                fields.append(f"{f}=None")
            elif isinstance(val, dict):
                fields.append(f"{f}=dict({len(val)} entries)")
            elif isinstance(val, (set, list)):
                fields.append(f"{f}={type(val).__name__}({len(val)} items)")
            else:
                fields.append(f"{f}={type(val).__name__}")
        methods = [
            name
            for name in dir(self)
            if not name.startswith("_")
            and callable(getattr(self, name))
            and name not in self.__dataclass_fields__
        ]
        fields.append(f"methods={methods}")
        body = ",\n  ".join(fields)
        return f"TaxonomicUtils(\n  {body}\n)"
        
    def load_a2t(self, accessions: List[str], low_memory: bool = None, extend: bool = False):
        """Load accession-to-taxon mappings, optionally extending the existing map."""
        if low_memory is None:
            low_memory = self._low_memory
        if extend:
            existing = dict(self.a2t or {})
            accessions = sorted({
                accession
                for accession in parse_accession(accessions, version=True)
                if accession not in existing
            })
            if not accessions:
                self.a2t = existing
                self._a2t_checked = True
                return
        self.a2t = build_a2t(
            accessions,
            low_memory=low_memory,
            verbose=not self._a2t_checked,
        )
        if extend:
            existing.update(self.a2t or {})
            self.a2t = existing
        self._a2t_checked = True

    def parse_accession(self, strings, version: bool = True):
        """Return accession IDs parsed from strings."""
        return parse_accession(strings, version=version)

    def get_rank_order(self):
        """Return canonical rank codes in taxonomic order."""
        return list(RANK_ORDER)

    def get_t2a(self, taxa, low_memory: bool = None):
        """Return accessions assigned to the provided taxa."""
        if low_memory is None:
            low_memory = self._low_memory
        accessions = get_t2a(
            taxa,
            low_memory=low_memory,
            verbose=not self._a2t_checked,
        )
        self._a2t_checked = True
        return accessions

    def _load_tree(self):
        tree = defaultdict(list)
        for k, v in self.parent.items():
            if v is not None:
                tree[int(v)].append(int(k))
        self._tree = tree

    def get_subtree(self, taxon):
        """Return all descendant taxa for a taxon, including the taxon itself."""
        if self._tree is None:
            self._load_tree()
        result = [taxon]
        if taxon in self._tree:
            for child in self._tree[taxon]:
                result.extend(get_subtree(child, self._tree))
        return result

    def sort_taxa(self, taxa):
        """Return unique taxa sorted in taxonomic order."""
        present = set(_as_taxa_list(taxa))
        rank = dict(zip(self.nodes["taxon"], self.nodes["rank_code"]))
        return taxonomic_order(present, self.parent, rank, self.names)

    def get_lca(self, a, b):
        """Return the lowest common ancestor taxon for two taxa."""
        return get_lca(int(a), int(b), self.parent)

    def get_branch(self, taxon):
        """Return the root-to-taxon branch for a taxon."""
        branch = []
        cur = int(taxon)
        seen = set()
        while cur is not None and cur not in seen:
            branch.append(cur)
            seen.add(cur)
            cur = self.parent.get(cur)
        return branch[::-1]

    def higher_than_rank(self, taxa, rank):
        """Return booleans indicating whether taxa are higher than the given rank."""
        rank_code = _rank_to_code(rank)
        threshold = RANK_ORDER[rank_code]
        rank_idx = dict(zip(self.nodes["taxon"], self.nodes["rank_idx"]))
        return np.array(
            [rank_idx.get(taxon, threshold) < threshold for taxon in _as_taxa_list(taxa)],
            dtype=bool,
        )

def download_taxonomy(
    accessions: List[str]=None,
    low_memory: bool=True,
    targets_json=None,
    rebuild: bool=False,
) -> TaxonomicUtils:
    """Download/load taxonomy resources and return a TaxonomicUtils object."""
    save_path = TAXUTILS_GLOBALS["save_folder"]
    os.makedirs(save_path, exist_ok=True)

    names_path = os.path.join(save_path, "names.dmp")
    nodes_path = os.path.join(save_path, "nodes.dmp")

    if rebuild or not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        logger.info(f"Downloading {names_path}, {nodes_path}...")
        tarball_path = os.path.join(save_path, "taxdump.tar.gz")
        url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
        _download_file(url, tarball_path)

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
    else:
        logger.info(
            "names.dmp and nodes.dmp exist in $TAXUTILS_GLOBALS "
            f"({save_path}), skipping download."
        )

    if targets_json is None:
        targets_json = os.path.join(save_path, "targets.json")
        if rebuild or not os.path.exists(targets_json):
            for url in TAXUTILS_GLOBALS["pathogen_dict_urls"]:
                try:
                    logger.info(f"Downloading targets.json from {url}...")
                    _download_file(url, targets_json)
                    break
                except Exception as e:
                    logger.warning(f"Failed to download from {url}: {e}")
            else:
                raise RuntimeError("Could not download targets.json from any URL.")
    logger.info(f"Building nodes...")
    names = build_names(names_path)
    nodes = build_nodes(nodes_path, names)
    parent = build_parent(nodes)
    target_taxa = build_target_taxa(nodes, names, targets_json=targets_json)
    if rebuild or not low_memory:
        _ensure_default_a2t_db(rebuild=rebuild)
    a2t = None
    if accessions is not None:
        a2t = build_a2t(accessions, low_memory=low_memory)
        a2t[TAXUTILS_GLOBALS["UNCLASSIFIED"]] = "unclassified"
    names[2697049] = "SARS-CoV-2"
    names[694009] = "SARS-related-CoV"
    return TaxonomicUtils(
        names=names,
        nodes=nodes,
        target_taxa=target_taxa,
        a2t=a2t,
        parent=parent,
        _low_memory=low_memory,
    )


def taxutils(
    accessions: List[str]=None,
    low_memory: bool=True,
    targets_json=None,
    rebuild: bool=False,
) -> TaxonomicUtils:
    """Build and return a TaxonomicUtils object."""
    return download_taxonomy(
        accessions=accessions,
        low_memory=low_memory,
        targets_json=targets_json,
        rebuild=rebuild,
    )

    
def _as_taxa_list(taxa):
    if np.isscalar(taxa):
        values = [taxa]
    elif isinstance(taxa, pd.Series):
        values = taxa.tolist()
    elif isinstance(taxa, np.ndarray):
        values = taxa.ravel().tolist()
    else:
        values = list(taxa)
    return [int(t) for t in values if not pd.isna(t)]


def _as_string_list(strings):
    if isinstance(strings, str):
        values = [strings]
    elif isinstance(strings, pd.Series):
        values = strings.tolist()
    elif isinstance(strings, np.ndarray):
        values = strings.ravel().tolist()
    else:
        values = list(strings)
    return [str(value) for value in values if not pd.isna(value)]


def parse_accession(strings, version: bool = True):
    """Return accession IDs parsed from strings."""
    accessions = []
    for text in _as_string_list(strings):
        for match in ACCESSION_PATTERN.finditer(text.upper()):
            accession = match.group(1)
            accession_version = match.group(2)
            if version and accession_version is not None:
                accession = f"{accession}.{accession_version}"
            accessions.append(accession)
    return accessions


def _rank_to_code(rank):
    key = str(rank).strip().upper()
    if key not in RANK_ALIASES:
        valid = ", ".join(RANK_ORDER)
        raise ValueError(f"rank must be one of: {valid}")
    return RANK_ALIASES[key]


def _ensure_a2t_db(gz_path, a2t_db, verbose=True, rebuild=False):
    """Build the SQLite a2t db from gz if it doesn't exist, ensuring both indexes exist."""
    if not os.path.exists(gz_path):
        gz_path = download_a2t(verbose=verbose)
    if rebuild and os.path.exists(a2t_db):
        os.remove(a2t_db)
    if not os.path.exists(a2t_db):
        logger.info("Building SQLite db from gz file, this will take a while...")
        conn = sqlite3.connect(a2t_db)
        cur = conn.cursor()
        cur.execute("CREATE TABLE a2t (accession TEXT, taxid INTEGER)")
        with gzip.open(gz_path, "rt") as f:
            next(f)
            batch = []
            for line in f:
                parts = line.strip().split("\t")
                batch.append((parts[1], int(parts[2])))
                if len(batch) == 100_000:
                    cur.executemany("INSERT INTO a2t VALUES (?, ?)", batch)
                    batch.clear()
            if batch:
                cur.executemany("INSERT INTO a2t VALUES (?, ?)", batch)
        cur.execute("CREATE INDEX idx_accession ON a2t (accession)")
        cur.execute("CREATE INDEX idx_taxid ON a2t (taxid)")
        conn.commit()
        conn.close()
        logger.info("SQLite db built.")
    else:
        # Ensure the taxon lookup index exists on dbs built before this index was added.
        conn = sqlite3.connect(a2t_db)
        cur = conn.cursor()
        cur.execute("CREATE INDEX IF NOT EXISTS idx_taxid ON a2t (taxid)")
        conn.commit()
        conn.close()


def _ensure_default_a2t_db(verbose=True, rebuild=False):
    gz_path = download_a2t(verbose=verbose, rebuild=rebuild)
    a2t_db = os.path.join(TAXUTILS_GLOBALS["save_folder"], "nucl_gb.accession2taxid.db")
    _ensure_a2t_db(gz_path, a2t_db, verbose=verbose, rebuild=rebuild)
    return gz_path, a2t_db

def download_a2t(verbose=True, rebuild=False):
    gz_path = os.path.join(TAXUTILS_GLOBALS["save_folder"], "nucl_gb.accession2taxid.gz")
    a2t_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
    if rebuild or not os.path.exists(gz_path):
        os.makedirs(os.path.dirname(gz_path), exist_ok=True)
        logger.info(f"Downloading {gz_path}...")
        _download_file(a2t_url, gz_path)
    else:
        if verbose:
            logger.info(f"{gz_path} already exists, skipping download.")

    return gz_path

def build_a2t(accessions, low_memory=True, verbose=True):
    accessions = parse_accession(accessions, version=True)

    if low_memory:
        gz_path = download_a2t(verbose=verbose)
        accession_set = set(accessions) if not isinstance(accessions, set) else accessions
        a2t = {}
        with gzip.open(gz_path, 'rt') as f:
            header = next(f).strip().split("\t")
            acc_idx = header.index("accession.version")
            taxon_idx = header.index("taxid")
            for line in f:
                parts = line.strip().split("\t")
                if parts[acc_idx] in accession_set:
                    a2t[parts[acc_idx]] = int(parts[taxon_idx])
                    if len(a2t) == len(accession_set):
                        break
        return a2t

    # SQLite path
    _, a2t_db = _ensure_default_a2t_db(verbose=verbose)

    conn = sqlite3.connect(a2t_db)
    acc_df = pd.DataFrame({"accession": list(accessions)})
    acc_df.to_sql("tmp_accs", conn, if_exists="replace", index=False)
    result = pd.read_sql("SELECT t.accession, a.taxid FROM tmp_accs t JOIN a2t a ON t.accession = a.accession", conn)
    conn.close()
    return dict(zip(result["accession"], result["taxid"]))


def get_t2a(taxa, low_memory=True, verbose=True):
    """Return the set of accessions belonging to the given taxa.

    When low_memory=True, scans the compressed file directly.
    When low_memory=False, uses a local SQLite database for speed.
    """
    if low_memory:
        gz_path = download_a2t(verbose=verbose)
        taxon_set = {int(t) for t in taxa}
        accessions = set()
        with gzip.open(gz_path, 'rt') as f:
            header = next(f).strip().split("\t")
            acc_idx = header.index("accession.version")
            taxon_idx = header.index("taxid")
            for line in f:
                parts = line.strip().split("\t")
                if int(parts[taxon_idx]) in taxon_set:
                    accessions.add(parts[acc_idx])
        return accessions

    # SQLite path
    _, a2t_db = _ensure_default_a2t_db(verbose=verbose)

    conn = sqlite3.connect(a2t_db)
    taxon_df = pd.DataFrame({"taxid": [int(t) for t in taxa]})
    taxon_df.to_sql("tmp_taxa", conn, if_exists="replace", index=False)
    result = pd.read_sql("SELECT accession FROM a2t JOIN tmp_taxa ON a2t.taxid = tmp_taxa.taxid", conn)
    conn.close()
    return set(result["accession"])

def taxonomic_order(present, parent, rank, names):
    anc, stack = set(), list(present)
    while stack:
        t = stack.pop()
        p = parent.get(t)
        if p is not None and p not in anc:
            anc.add(p)
            stack.append(p)

    nodes = present | anc
    children = {t: [] for t in nodes}
    for t in nodes:
        p = parent.get(t)
        if p in nodes:
            children[p].append(t)

    def child_key(t):
        return (str(rank.get(t, "")), str(names.get(t, "")), int(t))

    for k in children:
        children[k].sort(key=child_key)

    special_order = [0,1,9606,2,10239]
    roots = sorted(
        [t for t in nodes if parent.get(t) not in nodes],
        key=lambda t: (
            t not in special_order,
            special_order.index(t) if t in special_order else float("inf"),
            child_key(t),
        ),
    )

    order, seen = [], set()

    def dfs(u):
        if u in seen: return
        seen.add(u)
        if u in present:
            order.append(u)
        for v in children.get(u, []):
            dfs(v)

    for r in roots:
        dfs(r)
    for t in present:
        if t not in seen:
            order.append(t)
    return order


def build_parent(nodes):
    parent = dict(zip(nodes["taxon"], nodes["parent"]))
    parent[1] = None
    return parent

def build_names(names_path):
    names = {}
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxon = parts[0]
                name = parts[1]
                names[int(taxon)] = name
    names[TAXUTILS_GLOBALS["UNCLASSIFIED"]] = "unclassified" # -2
    return names


def get_subtree(taxon, tree):
    """Return all descendant taxa for a taxon from a parent-to-children tree."""
    result = [taxon]
    if taxon in tree:
        for child in tree[taxon]:
            result.extend(get_subtree(child, tree))
    return result


def assign_rank_codes(parent, rank_map):
    rank_code_cache = {}
    visiting = set()

    def code_base(code):
        return code[0]

    def code_depth(code):
        suffix = code[1:]
        return int(suffix) if suffix else 1

    def next_subrank(code):
        return f"{code_base(code)}{code_depth(code) + 1}"

    def rank_code(taxon):
        if taxon in rank_code_cache:
            return rank_code_cache[taxon]
        if taxon in visiting:
            rank_code_cache[taxon] = "R"
            return "R"

        visiting.add(taxon)
        if taxon == 0:
            code = "U"
        elif taxon == 1:
            code = "R"
        else:
            raw_code = MAJOR_RANK_TO_CODE.get(rank_map.get(taxon, ""))
            parent_taxon = parent.get(taxon)
            if parent_taxon is None or parent_taxon == taxon or parent_taxon not in parent:
                code = raw_code or "R"
            else:
                parent_code = rank_code(parent_taxon)
                parent_base = code_base(parent_code)
                if raw_code and RANK_ORDER[raw_code] > RANK_ORDER[parent_base]:
                    code = raw_code
                else:
                    code = next_subrank(parent_code)
        visiting.remove(taxon)
        rank_code_cache[taxon] = code
        return code

    return {taxon: rank_code(taxon) for taxon in parent}


def build_nodes(nodes_path, names):
    nodes = pd.read_csv(
        nodes_path, sep="|", header=None, usecols=[0,1,2],
        names=["taxon","parent","rank"], dtype={"taxon":int,"parent":int,"rank":str},
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

def get_parents(taxon, parent_map, rank_idx, rank="F"):
    parents = set()
    threshold = RANK_ORDER[_rank_to_code(rank)]
    cur_node = taxon
    while True:
        cur_node = parent_map.get(cur_node)
        if cur_node is None:
            break
        if rank_idx.get(cur_node, threshold - 1) < threshold:
            break
        parents.add(cur_node)
    return parents
    
def build_target_taxa(nodes, names, targets_json):
    with open(targets_json) as f:
        pdict = json.load(f)
    pathogen_taxa = {int(v) for v in pdict["pathogens"].values()}

    parent = build_parent(nodes)
    tree = defaultdict(list)
    for k, v in parent.items():
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            tree[int(v)].append(int(k))
    rank = dict(zip(nodes["taxon"], nodes["rank_code"]))
    rank_idx = dict(zip(nodes["taxon"], nodes["rank_idx"]))
    
    taxa = set()
    for taxon in pathogen_taxa:
        taxa.update(get_subtree(taxon, tree))
        taxa.update(get_parents(taxon, parent, rank_idx, rank="F"))

    return taxonomic_order(taxa, parent, rank, names)

def get_lca(a, b, parent_dict):
    if a == b:
        return a
    path_a = []
    cur = a
    while cur:
        path_a.append(cur)
        cur = parent_dict.get(cur)
    path_b = []
    cur = b
    while cur:
        path_b.append(cur)
        cur = parent_dict.get(cur)
    path_a = path_a[::-1]  # leaf to root reversed to root to leaf
    path_b = path_b[::-1]
    i = 0
    min_len = min(len(path_a), len(path_b))
    while i < min_len and path_a[i] == path_b[i]:
        i += 1
    return int(path_a[i-1]) if i > 0 else 1
