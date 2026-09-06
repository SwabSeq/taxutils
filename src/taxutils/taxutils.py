"""Taxonomy utilities backed by the native Rust extension.

`TaxonomicUtils` holds the loaded taxonomy and answers tree, rank and accession
questions. Resource loading lives in `resources`, tree algorithms in `tree`,
topology metrics in `topology`, and input coercion in `coerce`.

The save folder and the worker thread count are carried on the instance and
passed explicitly to every backend call; nothing reads process state at call
time.
"""

import os
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
import pandas as pd

from .coerce import as_taxa_list, is_taxon_scalar, pairwise_taxa_result, pairwise_values
from .parse import accessions_for_lookup, parse_accession
from .ranks import assign_rank_codes, rank_to_code
from .resources import (
    build_a2t,
    build_names,
    build_nodes,
    build_target_taxa,
    download_targets,
    download_taxdump,
    ensure_a2t_db,
    get_t2a,
)
from .topology import TopologyMixin
from .tree import (
    build_children,
    build_parent,
    get_lca,
    get_parents,
    get_subtree,
    taxonomic_order,
    tree_depth,
)
from .utils import (
    CANONICAL_RANK_NAMES,
    MAJOR_RANK_TO_CODE,
    RANK_ALIASES,
    RANK_ORDER,
    UNCLASSIFIED,
    get_logger,
    resolve_save_folder,
)

logger = get_logger(__name__)


@dataclass
class TaxonomicUtils(TopologyMixin):
    names: dict
    nodes: dict
    target_taxa: list
    a2t: dict = None
    parent: dict = None
    _low_memory: bool = True
    _wgs: bool = False
    _save_folder: Optional[str] = None
    _threads: Optional[int] = None

    def __post_init__(self):
        if self.parent is None:
            self.parent = build_parent(self.nodes)
        self.parent = {
            int(taxon): None if pd.isna(parent) else int(parent)
            for taxon, parent in self.parent.items()
        }
        self._save_folder = resolve_save_folder(self._save_folder)
        self._tree = None
        self._descendant_index = None
        self._depth = {}
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

    def load_a2t(
        self,
        accessions: List[str],
        low_memory: bool = None,
        extend: bool = False,
        wgs: bool = None,
    ):
        """Load accession-to-taxon mappings, optionally extending the existing map."""
        if low_memory is None:
            low_memory = self._low_memory
        if wgs is None:
            wgs = self._wgs
        if extend:
            existing = dict(self.a2t or {})
            accessions = sorted({
                accession
                for accession in accessions_for_lookup(accessions)
                if accession not in existing
            })
            if not accessions:
                self.a2t = existing
                self._a2t_checked = True
                return
        self.a2t = build_a2t(
            accessions,
            save_folder=self._save_folder,
            low_memory=low_memory,
            wgs=wgs,
            threads=self._threads,
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

    def get_t2a(self, taxa, low_memory: bool = None, wgs: bool = None):
        """Return accessions assigned directly to the provided taxa.

        By default, use the constructor's lookup mode. Pass low_memory=False
        to use the indexed SQLite database; low_memory=True scans the entire
        compressed accession mapping on every call, even for a single taxid.
        """
        if low_memory is None:
            low_memory = self._low_memory
        if wgs is None:
            wgs = self._wgs
        accessions = get_t2a(
            taxa,
            save_folder=self._save_folder,
            low_memory=low_memory,
            wgs=wgs,
            threads=self._threads,
        )
        self._a2t_checked = True
        return accessions

    def _load_tree(self):
        self._tree = build_children(self.parent)

    def _load_descendant_index(self):
        if self._tree is None:
            self._load_tree()

        start = {}
        end = {}
        visited = set()
        tick = 0
        roots = [
            taxon
            for taxon, parent in self.parent.items()
            if parent is None or parent == taxon or parent not in self.parent
        ]

        for root in roots + list(self.parent):
            if root in visited:
                continue
            stack = [(int(root), False)]
            while stack:
                node, exiting = stack.pop()
                if exiting:
                    end[node] = tick - 1
                    continue
                if node in visited:
                    continue
                visited.add(node)
                start[node] = tick
                tick += 1
                stack.append((node, True))
                for child in reversed(self._tree.get(node, [])):
                    if child not in visited:
                        stack.append((child, False))

        self._descendant_index = (start, end)

    def _get_depth(self, taxon):
        taxon = int(taxon)
        if taxon in self._depth:
            return self._depth[taxon]

        path = []
        seen = set()
        cur = taxon
        while cur is not None and cur not in seen and cur not in self._depth:
            path.append(cur)
            seen.add(cur)
            parent = self.parent.get(cur)
            cur = None if parent == cur else parent

        depth = self._depth.get(cur, -1)
        for node in reversed(path):
            depth += 1
            self._depth[node] = depth
        return self._depth[taxon]

    def get_subtree(self, taxon):
        """Return all descendant taxa for a taxon, including the taxon itself."""
        if self._tree is None:
            self._load_tree()
        result = [taxon]
        if taxon in self._tree:
            for child in self._tree[taxon]:
                result.extend(get_subtree(child, self._tree))
        return result

    def is_leaf(self, taxon):
        """Return whether taxa have no children in the taxonomy tree."""
        if self._tree is None:
            self._load_tree()

        def check(value):
            if pd.isna(value):
                return False
            return int(value) not in self._tree

        if isinstance(taxon, pd.Series):
            return taxon.map(check).astype(bool)

        if isinstance(taxon, np.ndarray):
            values = taxon.astype(object)
            return np.vectorize(check, otypes=[bool])(values)

        if np.isscalar(taxon):
            return check(taxon)

        return [check(value) for value in taxon]

    def is_child(self, taxon_a, taxon_b):
        """Return whether taxon_a is a direct child of taxon_b."""
        def check(a, b):
            if pd.isna(a) or pd.isna(b):
                return False
            return self.parent.get(int(a)) == int(b)

        return pairwise_taxa_result(taxon_a, taxon_b, check, name="is_child")

    def is_descendent(self, taxon_a, taxon_b):
        """Return whether taxon_a is a strict descendant of taxon_b."""
        if self._descendant_index is None:
            self._load_descendant_index()
        start, end = self._descendant_index

        def check(a, b):
            if pd.isna(a) or pd.isna(b):
                return False
            a = int(a)
            b = int(b)
            if a == b or a not in start or b not in start:
                return False
            return start[b] <= start[a] <= end[b]

        return pairwise_taxa_result(taxon_a, taxon_b, check, name="is_descendent")

    def _ancestor_at_rank(self, taxon, rank, rank_base=None):
        rank_code = rank_to_code(rank)
        if rank_base is None:
            rank_base = dict(zip(self.nodes["taxon"], self.nodes["rank_base"]))
        anchor = int(taxon)
        for branch_taxon in reversed(self.get_branch(taxon)):
            if rank_base.get(branch_taxon) == rank_code:
                return branch_taxon
        return anchor

    def get_ancestor(self, taxon, anchor_rank):
        """Return the nearest ancestor at a rank, preserving input type."""
        rank_code = rank_to_code(anchor_rank)
        rank_base = dict(zip(self.nodes["taxon"], self.nodes["rank_base"]))

        def get_one(value):
            if pd.isna(value):
                return np.nan
            return self._ancestor_at_rank(value, rank_code, rank_base)

        if isinstance(taxon, pd.Series):
            return taxon.map(get_one)

        if isinstance(taxon, np.ndarray):
            values = [get_one(value) for value in taxon.astype(object).ravel()]
            return np.asarray(values).reshape(taxon.shape)

        if np.isscalar(taxon):
            return get_one(taxon)

        return [get_one(value) for value in taxon]

    def sort_taxa(self, taxa):
        """Return unique taxa sorted in taxonomic order."""
        present = set(as_taxa_list(taxa))
        rank = dict(zip(self.nodes["taxon"], self.nodes["rank_code"]))
        return taxonomic_order(present, self.parent, rank, self.names)

    def format_tree(self, taxa, include_ancestors: bool = True, root: int = 1, indent: str = "\t"):
        """Return an indented taxonomic tree as a Series indexed by taxon."""
        taxa = set(as_taxa_list(taxa))
        tree_taxa = set()

        if include_ancestors:
            for taxon in taxa:
                cur = int(taxon)
                seen = set()
                while cur is not None and cur not in seen:
                    tree_taxa.add(cur)
                    if cur == root:
                        break
                    seen.add(cur)
                    cur = self.parent.get(cur)
        else:
            tree_taxa = set(taxa)

        order = self.sort_taxa(tree_taxa)
        depth = tree_depth(order, self.parent, root=root)
        names = {
            taxon: f"{indent * depth.get(taxon, 0)}{self.names.get(taxon, str(taxon))}"
            for taxon in order
        }
        return pd.Series(names, name="name").rename_axis("taxon")

    def get_lca(self, a, b):
        """Return the lowest common ancestor taxon for two taxa."""
        a = int(a)
        b = int(b)
        self._get_depth(a)
        self._get_depth(b)
        return get_lca(a, b, self.parent, self._depth)

    def get_distance(self, a, b):
        """Return edge distance between two taxa through their lowest common ancestor."""
        a = int(a)
        b = int(b)
        self._get_depth(a)
        self._get_depth(b)
        lca = get_lca(a, b, self.parent, self._depth)
        return self._depth.get(a, 0) + self._depth.get(b, 0) - 2 * self._depth.get(lca, 0)

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
        rank_code = rank_to_code(rank)
        threshold = RANK_ORDER[rank_code]
        rank_idx = dict(zip(self.nodes["taxon"], self.nodes["rank_idx"]))
        return np.array(
            [rank_idx.get(taxon, threshold) < threshold for taxon in as_taxa_list(taxa)],
            dtype=bool,
        )



def taxutils(
    accessions: List[str] = None,
    low_memory: bool = True,
    targets_json=None,
    wgs: bool = False,
    save_folder=None,
    threads: Optional[int] = None,
    refresh: bool = False,
) -> TaxonomicUtils:
    """Download or load taxonomy resources and return a `TaxonomicUtils`.

    Anything missing is downloaded and a missing or unusable accession
    database is built, so a first run needs no flags. `refresh=True` re-fetches
    the managed taxonomy files and brings an existing accession database up to
    date, applying only the rows NCBI added, changed or withdrew.

    `threads` is the worker count for every parallel stage, including the
    accession database build; `None` uses all logical CPUs. `save_folder`
    defaults to `$TAXUTILS_GLOBALS`, then `./taxutils/`.
    """
    save_path = resolve_save_folder(save_folder)
    os.makedirs(save_path, exist_ok=True)

    names_path = os.path.join(save_path, "names.dmp")
    nodes_path = os.path.join(save_path, "nodes.dmp")

    if refresh or not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        download_taxdump(save_path, names_path, nodes_path)
    else:
        logger.info(
            f"names.dmp and nodes.dmp exist in {save_path}, skipping download."
        )

    if targets_json is None:
        targets_json = os.path.join(save_path, "targets.json")
        if refresh or not os.path.exists(targets_json):
            download_targets(targets_json)

    logger.info("Building nodes...")
    names = build_names(names_path)
    nodes = build_nodes(nodes_path, names)
    parent = build_parent(nodes)
    target_taxa = build_target_taxa(nodes, names, targets_json=targets_json)

    if refresh or not low_memory:
        ensure_a2t_db(
            save_folder=save_path,
            wgs=wgs,
            refresh=refresh,
            threads=threads,
        )

    a2t = None
    if accessions is not None:
        a2t = build_a2t(
            accessions,
            save_folder=save_path,
            low_memory=low_memory,
            wgs=wgs,
            threads=threads,
        )
        a2t[UNCLASSIFIED] = "unclassified"

    names[2697049] = "SARS-CoV-2"
    names[694009] = "SARS-related-CoV"
    return TaxonomicUtils(
        names=names,
        nodes=nodes,
        target_taxa=target_taxa,
        a2t=a2t,
        parent=parent,
        _low_memory=low_memory,
        _wgs=wgs,
        _save_folder=save_path,
        _threads=threads,
    )


# `download_taxonomy` was the original entry point and stays as an alias.
download_taxonomy = taxutils
