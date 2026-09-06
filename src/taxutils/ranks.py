"""Rank code normalization and NCBI rank correction.

NCBI ranks are irregular: many nodes carry "no rank" or "clade", and a child can
repeat its parent's rank. `assign_rank_codes` walks each lineage and assigns a
canonical code, appending a depth suffix (``G2``) whenever a node cannot claim a
rank strictly below its parent's.
"""

from .utils import MAJOR_RANK_TO_CODE, RANK_ALIASES, RANK_ORDER


def rank_to_code(rank):
    """Return the canonical single-letter code for a rank name or code."""
    key = str(rank).strip().upper()
    if key not in RANK_ALIASES:
        valid = ", ".join(RANK_ORDER)
        raise ValueError(f"rank must be one of: {valid}")
    return RANK_ALIASES[key]


def assign_rank_codes(parent, rank_map):
    """Return a corrected rank code for every taxon in `parent`."""
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
