# taxutils.py

import pandas as pd
import numpy as np
from collections import defaultdict
from typing import List
from dataclasses import dataclass
import os, json, urllib.request, tarfile, gzip, re
import sqlite3, subprocess

try:
    from .utils import TAXUTILS_GLOBALS, get_logger
except ImportError:
    from utils import TAXUTILS_GLOBALS, get_logger

logger = get_logger(__name__)

@dataclass
class TaxonomicUtils:
    names: dict
    nodes: dict
    target_taxids: set
    a2t: dict = None
    tree: defaultdict = None
    _is_matched: bool = False
    parent: dict = None
    
    def __repr__(self):
        fields = []
        for f in self.__dataclass_fields__:
            val = getattr(self, f)
            if val is None:
                fields.append(f"{f}=None")
            elif isinstance(val, dict):
                fields.append(f"{f}=dict({len(val)} entries)")
            elif isinstance(val, (set, list)):
                fields.append(f"{f}={type(val).__name__}({len(val)} items)")
            else:
                fields.append(f"{f}={type(val).__name__}")
        body = ",\n  ".join(fields)
        return f"TaxonomicUtils(\n  {body}\n)"
        
    def load_a2t(self, accessions: List[str], sqlite: bool = True):
        self.a2t = build_a2t(accessions, sqlite=sqlite)

    def load_tree(self):
        parent = dict(zip(
            self.nodes["taxid"],
            self.nodes["parent"])
        )
        parent[1] = None
        tree = defaultdict(list)
        for k, v in parent.items():
            if v is not None:
                tree[int(v)].append(int(k))
        self.tree = tree

    def get_subtree(self, taxid):
        if self.tree is None:
            self.load_tree()
        result = [taxid]
        if taxid in self.tree:
            for child in self.tree[taxid]:
                result.extend(get_subtree(child, self.tree))
        return result

    def match_library(self, config):
        logger.info("Matching library to target taxa.")
        acc_ids = extract_accession_ids(config.library)
        self.load_a2t(acc_ids, sqlite=config.sqlite)
        library_taxids = set()
        for acc_id in acc_ids:
            add_to_lib = self.a2t.get(acc_id, None)
            if add_to_lib is not None:
                library_taxids.add(add_to_lib)
        
        if self.tree is None:
            self.load_tree()
        parent = dict(zip(self.nodes['taxid'], self.nodes['parent']))
        parent[1] = None
        higher_than_F = dict(zip(self.nodes['taxid'], self.nodes['higher_than_F']))
        
        target_set = set(self.target_taxids)
        matched = target_set.intersection(library_taxids)
        
        for tid in library_taxids:
            cur = parent.get(tid)
            while cur is not None and not higher_than_F.get(cur, True):
                if cur in target_set:
                    matched.add(cur)
                cur = parent.get(cur)
        
        self.target_taxids = matched
        self._is_matched = True

    def load_parent(self):
        self.parent = dict(zip(self.nodes["taxid"], self.nodes["parent"]))

def extract_accession_ids(fasta_path):
    """Extract accession IDs with version numbers from FASTA headers using grep."""
    result = subprocess.run(
        ['grep', '-o', '^>.*', fasta_path],
        capture_output=True, text=True
    )
    
    headers = result.stdout
    accession_ids = re.findall(r'[A-Z]{1,2}_?\d+\.\d+', headers)
    
    return accession_ids

def download_taxonomy(accessions: List[str]=None, sqlite: bool=True, pathogen_json=None) -> TaxonomicUtils:
    save_path = TAXUTILS_GLOBALS["save_folder"]
    os.makedirs(save_path, exist_ok=True)

    names_path = os.path.join(save_path, "names.dmp")
    nodes_path = os.path.join(save_path, "nodes.dmp")

    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        logger.info(f"Downloading {names_path}, {nodes_path}...")
        tarball_path = os.path.join(save_path, "taxdump.tar.gz")
        url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
        urllib.request.urlretrieve(url, tarball_path)

        with tarfile.open(tarball_path, "r:gz") as tar:
            for member in tar.getmembers():
                if member.name in ["names.dmp", "nodes.dmp"]:
                    tar.extract(member, path=save_path)

        os.remove(tarball_path)
    else:
        logger.info(f"names.dmp and nodes.dmp exist in {save_path}, skipping download.")

    if pathogen_json is None:
        pathogen_json = os.path.join(save_path, "pathogen_dict.json")
        if not os.path.exists(pathogen_json):
            for url in TAXUTILS_GLOBALS["pathogen_dict_urls"]:
                try:
                    logger.info(f"Downloading pathogen_dict.json from {url}...")
                    urllib.request.urlretrieve(url, pathogen_json)
                    break
                except Exception as e:
                    logger.warning(f"Failed to download from {url}: {e}")
            else:
                raise RuntimeError("Could not download pathogen_dict.json from any URL.")
    logger.info(f"Building nodes...")
    names = build_names(names_path)
    nodes = build_nodes(nodes_path, names)
    target_taxids = build_target_taxids(
        nodes, names, pathogen_json=pathogen_json, extra_taxids=(9606,)
    )
    a2t = None
    if accessions is not None:
        a2t = build_a2t(accessions, sqlite=sqlite)
        a2t[TAXUTILS_GLOBALS["UNCLASSIFIED"]] = "unclassified"
        a2t[TAXUTILS_GLOBALS["UNMAPPED"]] = "unmapped"
    names[2697049] = "SARS-CoV-2"
    names[694009] = "SARS-related-CoV"
    return TaxonomicUtils(names=names, nodes=nodes, target_taxids=target_taxids, a2t=a2t)


def taxutils(accessions: List[str]=None, sqlite: bool=True, pathogen_json=None) -> TaxonomicUtils:
    return download_taxonomy(accessions=accessions, sqlite=sqlite, pathogen_json=pathogen_json)


TaxonomicData = TaxonomicUtils

    
def _ensure_a2t_db(gz_path, a2t_db):
    """Build the SQLite a2t db from gz if it doesn't exist, ensuring both indexes exist."""
    if not os.path.exists(gz_path):
        gz_path = download_a2t()
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
        # Ensure taxid index exists on dbs built before this index was added
        conn = sqlite3.connect(a2t_db)
        cur = conn.cursor()
        cur.execute("CREATE INDEX IF NOT EXISTS idx_taxid ON a2t (taxid)")
        conn.commit()
        conn.close()

def download_a2t():
    gz_path = os.path.join(TAXUTILS_GLOBALS["save_folder"], "nucl_gb.accession2taxid.gz")
    a2t_url = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"
    if not os.path.exists(gz_path):
        os.makedirs(os.path.dirname(gz_path), exist_ok=True)
        logger.info(f"Downloading {gz_path}...")
        urllib.request.urlretrieve(a2t_url, gz_path)
    else:
        logger.info(f"{gz_path} already exists, skipping download.")

    return gz_path

def build_a2t(accessions, sqlite=True):
    gz_path = download_a2t()

    if not sqlite:
        accession_set = set(accessions) if not isinstance(accessions, set) else accessions
        a2t = {}
        with gzip.open(gz_path, 'rt') as f:
            header = next(f).strip().split("\t")
            acc_idx = header.index("accession.version")
            taxid_idx = header.index("taxid")
            for line in f:
                parts = line.strip().split("\t")
                if parts[acc_idx] in accession_set:
                    a2t[parts[acc_idx]] = int(parts[taxid_idx])
                    if len(a2t) == len(accession_set):
                        break
        return a2t

    # sqlite path
    a2t_db = os.path.join(TAXUTILS_GLOBALS["save_folder"], "nucl_gb.accession2taxid.db")
    _ensure_a2t_db(gz_path, a2t_db)

    conn = sqlite3.connect(a2t_db)
    acc_df = pd.DataFrame({"accession": list(accessions)})
    acc_df.to_sql("tmp_accs", conn, if_exists="replace", index=False)
    result = pd.read_sql("SELECT t.accession, a.taxid FROM tmp_accs t JOIN a2t a ON t.accession = a.accession", conn)
    conn.close()
    return dict(zip(result["accession"], result["taxid"]))


def get_accessions_for_taxids(taxids, sqlite=True):
    """Return the set of accessions belonging to the given taxids.

    When sqlite=True (default), queries the local SQLite db for speed.
    When sqlite=False, scans the compressed gz file (keeps the file compressed
    but is much slower).
    """
    gz_path = download_a2t()

    if not sqlite:
        taxid_set = {int(t) for t in taxids}
        accessions = set()
        with gzip.open(gz_path, 'rt') as f:
            header = next(f).strip().split("\t")
            acc_idx = header.index("accession.version")
            taxid_idx = header.index("taxid")
            for line in f:
                parts = line.strip().split("\t")
                if int(parts[taxid_idx]) in taxid_set:
                    accessions.add(parts[acc_idx])
        return accessions

    # sqlite path
    a2t_db = os.path.join(TAXUTILS_GLOBALS["save_folder"], "nucl_gb.accession2taxid.db")
    _ensure_a2t_db(gz_path, a2t_db)

    taxid_list = [int(t) for t in taxids]
    conn = sqlite3.connect(a2t_db)
    taxid_df = pd.DataFrame({"taxid": [int(t) for t in taxids]})
    taxid_df.to_sql("tmp_taxids", conn, if_exists="replace", index=False)
    result = pd.read_sql("SELECT accession FROM a2t JOIN tmp_taxids ON a2t.taxid = tmp_taxids.taxid", conn)
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

def get_parent_tree(taxonomic_data):
    parent = dict(zip(
        taxonomic_data.nodes["taxid"],
        taxonomic_data.nodes["parent"])
    )
    parent[1] = None
    tree = defaultdict(list)
    for k, v in parent.items():
        if v is not None:
            tree[int(v)].append(int(k))
    return parent, tree

def build_names(names_path):
    names = {}
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 4 and parts[3] == "scientific name":
                taxid = parts[0]
                name = parts[1]
                names[int(taxid)] = name
    names[TAXUTILS_GLOBALS["UNMAPPED"]] = "unmapped" # -1
    names[TAXUTILS_GLOBALS["UNCLASSIFIED"]] = "unclassified" # -2
    return names

def get_subtree(taxid, tree):
    """
    Get all descendant taxids including itself.
    """
    result = [taxid]
    if taxid in tree:
        for child in tree[taxid]:
            result.extend(get_subtree(child, tree))
    return result

def rank_below(r):
    order = [
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
    ]
    rank_index = {r: i for i, r in enumerate(order)}
    if r == "species": return "species"
    if r not in rank_index:
        return None
    i = rank_index[r]
    if i+1 < len(order):
        return order[i+1]
    return None
    
def correct_rank(taxid, rank_map, parent, seen):
    if taxid in seen:
        return "root"
    seen.add(taxid)
    cur_rank = rank_map.get(taxid)
    if cur_rank is None or cur_rank == "root":
        return "root"
    if cur_rank != "no rank":
        return cur_rank
    return correct_rank(parent[taxid], rank_map, parent, seen)

def build_nodes(nodes_path, names):
    major_letters = ['U','R','D','K','P','C','O','F','G','S']
    major_order   = {'U':0,'R':1,'D':2,'K':3,'P':4,'C':5,'O':6,'F':7,'G':8,'S':9}
    major_rank_to_code = {
        "root":"R",
        "acellular root":"R",
        "cellular root":"R",
        "no rank":"NR",
        "clade":"C",
        "subfamily":"F",
        "domain":"D",
        "realm":"D",
        "kingdom":"K",
        "phylum":"P",
        "class":"C",
        "order":"O",
        "family":"F",
        "genus":"G",
        "species":"S",
    }
    
    nodes = pd.read_csv(
        nodes_path, sep="|", header=None, usecols=[0,1,2],
        names=["taxid","parent","rank"], dtype={"taxid":int,"parent":int,"rank":str},
        engine="python"
    )
    nodes["rank"] = nodes["rank"].str.strip().str.lower()
    # nodes["name"] = nodes["taxid"].map(names)
    nodes.loc[nodes["taxid"]==1, "parent"]
    parent = dict(zip(nodes["taxid"], nodes["parent"]))
    rank_map = dict(zip(nodes["taxid"], nodes["rank"]))
    
    _rank_code_cache = {}
    def _rank_code(t):
        if t in _rank_code_cache: return _rank_code_cache[t]
        if t == 0: _rank_code_cache[t] = "U"; return "U"
        if t == 1: _rank_code_cache[t] = "R"; return "R"
        steps, cur = 0, t
        while True:
            r = rank_map.get(cur, "")
            b = major_rank_to_code.get(r)
            if b:
                c = b if steps == 0 else f"{b}{steps}"
                _rank_code_cache[t] = c
                return c
            nxt = parent.get(cur)
            if nxt is None or nxt == cur: _rank_code_cache[t] = None; return None
            cur = nxt
            steps += 1

    nodes["rank"] = nodes["taxid"].apply(
        lambda t: correct_rank(t, rank_map, parent, set())
    )
    nodes["rank_code"] = nodes["taxid"].apply(_rank_code)
    nodes["rank_code"] = nodes["rank_code"]
    nodes["rank_base"] = nodes["rank_code"].str[0]
    nodes["rank_idx"] = nodes["rank_base"].map(major_order)

    for L in major_letters:
        ti = major_order[L]
        nodes[f"higher_than_{L}"] = nodes["rank_idx"] <  ti
    
    return nodes

def get_parents(tid, parent_map, higher_than_F):
    parents = set()
    cur_node = tid
    while True:
        cur_node = parent_map.get(cur_node)
        if cur_node is None:
            break
        if higher_than_F.get(cur_node, True):
            break
        parents.add(cur_node)
    return parents
    
def build_target_taxids(nodes, names, pathogen_json, extra_taxids=(9606,)):
    with open(pathogen_json) as f:
        pdict = json.load(f)
    pathogen_taxids = {int(v) for v in pdict["pathogens"].values()}

    parent = dict(zip(nodes["taxid"], nodes["parent"]))
    parent[1] = None
    tree = defaultdict(list)
    for k, v in parent.items():
        if v is not None and not (isinstance(v, float) and np.isnan(v)):
            tree[int(v)].append(int(k))
    rank = dict(zip(nodes["taxid"], nodes["rank"]))
    higher_than_F = dict(zip(nodes["taxid"], nodes["higher_than_F"]))
    
    taxids = set()
    for tid in pathogen_taxids:
        taxids.update(get_subtree(tid, tree))
        taxids.update(get_parents(tid, parent, higher_than_F))

    if extra_taxids:
        taxids.update(extra_taxids)
    return taxonomic_order(taxids, parent, rank, names)

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
