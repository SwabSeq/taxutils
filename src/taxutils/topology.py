"""Subtree topology metrics.

These are separated from `TaxonomicUtils` only to keep that class readable; the
mixin is applied to it and the methods behave exactly as before. `_topology_context`
memoizes per-anchor state so a batch of taxa sharing an anchor walks the subtree once.
"""

import numpy as np
import pandas as pd

from .coerce import as_taxa_list as _as_taxa_list
from .ranks import rank_to_code as _rank_to_code


class TopologyMixin:
    def topology(self, taxon, anchor_rank=None, stat=None):
        """Return subtree topology metrics or a single topology statistic."""
        stat = None if stat in (None, "") else str(stat)
        rank_code_map = (
            dict(zip(self.nodes["taxon"], self.nodes["rank_code"]))
            if stat is None
            else None
        )
        rank_base = (
            dict(zip(self.nodes["taxon"], self.nodes["rank_base"]))
            if anchor_rank is not None
            else None
        )
        needs_subtree_set = stat is None or stat in {
            "n_leaves",
            "max_children",
            "branching_taxa_fraction",
            "top_child_fraction",
        }

        if not np.isscalar(taxon):
            taxa = _as_taxa_list(taxon)
            context_cache = {} if anchor_rank is not None else None
            if stat is not None:
                values = []
                for value in taxa:
                    context = self._topology_context(
                        value,
                        anchor_rank=anchor_rank,
                        rank_base=rank_base,
                        needs_subtree_set=needs_subtree_set,
                        context_cache=context_cache,
                    )
                    values.append(self._topology_stat_value(context, stat))
                return pd.Series(values, index=taxa, name=stat).rename_axis("taxon")

            rows = []
            for value in taxa:
                context = self._topology_context(
                    value,
                    anchor_rank=anchor_rank,
                    rank_base=rank_base,
                    needs_subtree_set=needs_subtree_set,
                    context_cache=context_cache,
                )
                rows.append(self._topology_profile_from_context(
                    context,
                    rank_code_map=rank_code_map,
                ))
            return pd.DataFrame(rows)

        context = self._topology_context(
            taxon,
            anchor_rank=anchor_rank,
            rank_base=rank_base,
            needs_subtree_set=needs_subtree_set,
        )
        if stat is not None:
            return self._topology_stat_value(context, stat)
        return self._topology_profile_from_context(context, rank_code_map=rank_code_map)

    def _topology_stat_value(self, context, stat):
        stats = {
            "n_taxa": self._topology_n_taxa,
            "n_leaves": self._topology_n_leaves,
            "max_depth": self._topology_max_depth,
            "mean_depth": self._topology_mean_depth,
            "topology_scale": self._topology_scale,
            "max_children": self._topology_max_children,
            "branching_taxa_fraction": self._topology_branching_taxa_fraction,
            "top_child_fraction": self._topology_top_child_fraction,
        }
        if stat not in stats:
            valid = ", ".join(stats)
            raise ValueError(f"stat must be one of: {valid}")
        return stats[stat](context)

    def _topology_context(
        self,
        taxon,
        anchor_rank=None,
        rank_base=None,
        needs_subtree_set=True,
        context_cache=None,
    ):
        taxon = int(taxon)
        rank_code = None if anchor_rank is None else _rank_to_code(anchor_rank)

        if self._tree is None:
            self._load_tree()

        anchor = taxon if rank_code is None else self._ancestor_at_rank(
            taxon,
            rank_code,
            rank_base,
        )
        if context_cache is not None and anchor in context_cache:
            context = context_cache[anchor]
            context["taxon"] = taxon
            return context

        subtree = [int(node) for node in self.get_subtree(anchor)]
        context = {
            "taxon": taxon,
            "anchor": anchor,
            "subtree": subtree,
        }
        if needs_subtree_set:
            context["subtree_set"] = set(subtree)
        if context_cache is not None:
            context_cache[anchor] = context
        return context

    def _topology_relative_depths(self, context):
        if "relative_depths" in context:
            return context["relative_depths"]
        anchor_depth = self._get_depth(context["anchor"])
        context["relative_depths"] = [
            max(self._get_depth(node) - anchor_depth, 0)
            for node in context["subtree"]
        ]
        return context["relative_depths"]

    def _topology_child_counts(self, context):
        if "child_counts" in context:
            return context["child_counts"]
        subtree_set = context.setdefault("subtree_set", set(context["subtree"]))
        context["child_counts"] = [
            sum(child in subtree_set for child in self._tree.get(node, []))
            for node in context["subtree"]
        ]
        return context["child_counts"]

    def _topology_n_taxa(self, context):
        return len(context["subtree"])

    def _topology_n_leaves(self, context):
        return sum(count == 0 for count in self._topology_child_counts(context))

    def _topology_max_depth(self, context):
        relative_depths = self._topology_relative_depths(context)
        return max(relative_depths) if relative_depths else 0

    def _topology_mean_depth(self, context):
        relative_depths = self._topology_relative_depths(context)
        return float(np.mean(relative_depths)) if relative_depths else 0

    def _topology_scale(self, context):
        relative_depths = self._topology_relative_depths(context)
        descendant_depths = sorted(depth for depth in relative_depths if depth > 0)
        p95_depth = (
            descendant_depths[int(0.95 * (len(descendant_depths) - 1))]
            if descendant_depths
            else 0
        )
        return max(p95_depth, 1)

    def _topology_max_children(self, context):
        child_counts = self._topology_child_counts(context)
        return max(child_counts) if child_counts else 0

    def _topology_branching_taxa_fraction(self, context):
        n_taxa = self._topology_n_taxa(context)
        if not n_taxa:
            return 0
        child_counts = self._topology_child_counts(context)
        return sum(count > 0 for count in child_counts) / n_taxa

    def _topology_top_child_fraction(self, context):
        subtree_set = context["subtree_set"]
        subtree_sizes = {}
        for node in reversed(context["subtree"]):
            subtree_sizes[node] = 1 + sum(
                subtree_sizes[child]
                for child in self._tree.get(node, [])
                if child in subtree_set
            )

        immediate_children = [
            child for child in self._tree.get(context["anchor"], []) if child in subtree_set
        ]
        immediate_child_sizes = [subtree_sizes[child] for child in immediate_children]
        total_child_size = sum(immediate_child_sizes)
        return max(immediate_child_sizes) / total_child_size if total_child_size else 1

    def _topology_profile_from_context(self, context, rank_code_map):
        taxon = context["taxon"]
        anchor = context["anchor"]
        profile = pd.Series({
            "taxon": taxon,
            "name": self.names.get(taxon, str(taxon)),
            "rank_code": rank_code_map.get(taxon),
            "anchor_taxon": anchor,
            "anchor_name": self.names.get(anchor, str(anchor)),
            "anchor_rank_code": rank_code_map.get(anchor),
            "n_taxa": self._topology_n_taxa(context),
            "n_leaves": self._topology_n_leaves(context),
            "max_depth": self._topology_max_depth(context),
            "mean_depth": self._topology_mean_depth(context),
            "topology_scale": self._topology_scale(context),
            "max_children": self._topology_max_children(context),
            "branching_taxa_fraction": self._topology_branching_taxa_fraction(context),
            "top_child_fraction": self._topology_top_child_fraction(context),
        }, dtype=object)
        return profile

