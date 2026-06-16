# taxutils

Utilities for working with NCBI taxonomic data, accession-to-taxon mappings, taxonomy branches, corrected ranks, and pathogen target taxa.

# Install

### Conda Formula
```bash
conda install bioconda::taxutils
```

### Pip Formula
```bash
pip install taxutils
```

# Setup

`taxutils` stores downloaded taxonomy files (`names.dmp`, `nodes.dmp`), pathogen target metadata, and accession-to-taxon mappings in a global save directory. Set `TAXUTILS_GLOBALS` before importing the package if you want to control where these files live:

```bash
export TAXUTILS_GLOBALS=/path/to/taxutils/saves
```

If `TAXUTILS_GLOBALS` is not set, `taxutils` defaults to `./taxutils/` in the current working directory.

The first run downloads NCBI taxonomy files. Accession lookups also use the NCBI accession-to-taxon mapping, which is large. By default, `taxutils` uses low-memory mode and scans the compressed mapping directly. For faster repeated lookups, use `low_memory=False` to build or reuse a local SQLite database:

```python
from taxutils import taxutils

tu = taxutils(low_memory=False)
```

# Core usage

Core functions are listed here. See the example notebook for a fuller walkthrough.

```python
# Build object
tu = taxutils(accessions=None, low_memory=True, targets_json=None, rebuild=False)  # Build the taxonomy utility object.

# Accession parsing and mapping
tu.parse_accession(header_strings, version=True)          # Extract one accession per string.
tu.load_a2t(accessions, low_memory=None, extend=False)    # Load accession-to-taxon mappings into tu.a2t.
tu.get_t2a(taxa, low_memory=None)                         # Return accessions assigned to taxa.

# Tree queries
tu.get_branch(taxon)                                      # Return the root-to-taxon branch.
tu.get_subtree(taxon)                                     # Return taxon plus all descendants.
tu.get_lca(taxon_a, taxon_b)                              # Return the lowest common ancestor.
tu.get_distance(taxon_a, taxon_b)                         # Return tree edge distance through the LCA.
tu.sort_taxa(taxa)                                        # Sort taxa in taxonomic order.
tu.format_tree(taxa)                                      # Return an indented tree Series.
tu.topology(taxon, anchor_rank=None)                      # Return subtree topology metrics.
tu.topology(taxon, anchor_rank=None, stat="topology_scale")  # Return one topology statistic.

# Rank utilities
tu.get_rank_order()                                       # Return canonical rank codes.
tu.higher_than_rank(taxa, rank)                           # Test whether taxa are higher than a rank.
```

In taxutils, `accessions=list/of/accessions` can be passed to call load_a2t on construction of the taxutils object. A custom targets_json can similarly be passed in lieu of the default json explained below. `rebuild=True` redownloads the managed taxonomy, target, and accession files and rebuilds the SQLite accession database. `load_a2t` overwrites `tu.a2t` by default; pass `extend=True` to add missing mappings without discarding existing ones. Method-level `low_memory=None` uses the mode set when `tu` was built.

`parse_accession` accepts strings, lists, arrays, and pandas Series. It returns the first accession found from each string using the same container type where possible; missing accessions are returned as `"NA"`.

`get_lca(a, b)` returns the lowest common ancestor of two taxa. `get_distance(a, b)` returns the edge distance between two taxa through their lowest common ancestor. Depths are cached lazily as these methods are called.

`topology(taxon, anchor_rank=None, stat=None)` returns subtree topology metrics such as taxon count, leaf fraction, depth, branchiness, and `topology_scale`. Pass `anchor_rank="F"` to summarize the nearest family-level ancestor. With `stat=None`, a single taxon returns a Series and a list, array, or Series returns a DataFrame.

Pass `stat` to return one topology metric. A single taxon returns a scalar; a list, array, or Series returns a Series indexed by taxon.

## Topology columns

`tu.topology(...)` returns these columns:

- `taxon`: input taxon.
- `name`: input taxon name.
- `rank_code`: corrected rank code for the input taxon.
- `anchor_taxon`: subtree root used for the topology summary.
- `anchor_name`: anchor taxon name.
- `anchor_rank_code`: corrected rank code for the anchor.
- `n_taxa`: number of taxa in the anchor subtree.
- `n_leaves`: number of terminal taxa in the anchor subtree.
- `max_depth`: maximum number of edges below the anchor.
- `mean_depth`: average number of edges below the anchor.
- `topology_scale`: 95th percentile descendant depth, with minimum value 1.
- `max_children`: largest number of direct children from any node in the subtree.
- `branching_taxa_fraction`: fraction of subtree taxa with at least one child.
- `top_child_fraction`: fraction of the anchor subtree contained in its largest immediate child branch.

`tu.target_taxa` contains the default pathogen-derived target taxa. Use it directly for target filtering or movement checks.

# Rank correction

`taxutils` keeps the raw NCBI rank in `rank` and adds corrected rank columns. Canonical ranks (`R`, `D`, `K`, `P`, `C`, `O`, `F`, `G`, `S`) are used as anchors only when they move deeper than the corrected parent rank. Noncanonical ranks such as `no rank`, `clade`, and other unusual labels inherit position from the tree. If a child would be ranked at the same or a higher level than its parent, it is assigned a subrank such as `S2`, `S3`, or `F2`. The canonical name for the corrected rank is stored in `new_rank`.

# Target taxa

In ZarLab, we are working on metagenomics in the clinical setting, with the goal of creating an "agnostic diagnostic". We often want to look at broad array of taxa (`tu.target_taxa`) that could cause harm to people. In June 2024, CZI did the work of compiling a list of pathogenic taxa. I did the easy work of turning this into a json and uploading it to my website, so that it is available and easily accessed for all time (in case that link ever breaks). taxutils will extend the taxa list to include subtrees of each of those pathogenic taxa. It will additionally include SARS-CoV2, since it was excluded from CZI's list. If you find any other obvious, missing pathogens, please send me a note, so I can update my json. You can also update the target_taxa member variable yourself, or store an entirely different set of targets, if you wanted.

# Contact

Author: Will O'Brien  
Affiliation: Computer Science Department, UCLA  
Email: wob@cs.ucla.edu
