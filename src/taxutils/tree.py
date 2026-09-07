"""Pure taxonomy-tree algorithms.

Nothing here touches the filesystem, the network or the Rust backend: these
operate on plain parent/children dictionaries so they can be tested in
isolation and reused by both the class methods and the module-level functions.
"""

from collections import defaultdict

import numpy as np

from .utils import RANK_ORDER


def build_parent(nodes):
    parent = dict(zip(nodes["taxon"], nodes["parent"]))
    parent[1] = None
    return parent


def build_children(parent):
    """Return a parent-to-children map, skipping missing parents."""
    tree = defaultdict(list)
    for taxon, parent_taxon in parent.items():
        if parent_taxon is None:
            continue
        if isinstance(parent_taxon, float) and np.isnan(parent_taxon):
            continue
        tree[int(parent_taxon)].append(int(taxon))
    return tree


def get_subtree(taxon, tree):
    """Return all descendant taxa for a taxon from a parent-to-children tree."""
    result = []
    stack = [taxon]
    while stack:
        node = stack.pop()
        result.append(node)
        children = tree.get(node, [])
        stack.extend(reversed(children))
    return result


def get_parents(taxon, parent_map, rank_idx, rank="F"):
    from .ranks import rank_to_code

    parents = set()
    threshold = RANK_ORDER[rank_to_code(rank)]
    cur_node = taxon
    while True:
        cur_node = parent_map.get(cur_node)
        if cur_node is None:
            break
        if rank_idx.get(cur_node, threshold - 1) < threshold:
            break
        parents.add(cur_node)
    return parents


def get_lca(a, b, parent_dict, depth=None):
    """Return the lowest common ancestor of two taxa.

    With `depth` supplied the two lineages are levelled first and walked in
    step; without it both root paths are built and compared from the root.
    """
    if a == b:
        return a
    if depth is None:
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
        path_a = path_a[::-1]
        path_b = path_b[::-1]
        i = 0
        min_len = min(len(path_a), len(path_b))
        while i < min_len and path_a[i] == path_b[i]:
            i += 1
        return int(path_a[i - 1]) if i > 0 else 1

    a = int(a)
    b = int(b)
    depth_a = depth.get(a, 0)
    depth_b = depth.get(b, 0)

    while depth_a > depth_b:
        a = parent_dict.get(a)
        depth_a -= 1
    while depth_b > depth_a:
        b = parent_dict.get(b)
        depth_b -= 1

    while a != b:
        a = parent_dict.get(a)
        b = parent_dict.get(b)
        if a is None or b is None:
            return 1
    return int(a)


def tree_depth(order, parent, root=1):
    """Return display depth per taxon, counting only taxa present in `order`."""
    visible = set(order)
    depth = {}
    for taxon in order:
        chain = []
        cur = int(taxon)
        seen = set()
        while cur is not None and cur not in seen:
            chain.append(cur)
            if cur == root:
                break
            seen.add(cur)
            cur = parent.get(cur)

        visible_ancestors = [node for node in chain[::-1] if node in visible]
        for idx, node in enumerate(visible_ancestors):
            depth.setdefault(node, idx)

    return depth


def taxonomic_order(present, parent, rank, names):
    """Return `present` ordered by a depth-first walk of the taxonomy."""
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

    special_order = [0, 1, 9606, 2, 10239]
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
        if u in seen:
            return
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
