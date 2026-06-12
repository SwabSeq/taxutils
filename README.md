# taxutils

Utilities for working with NCBI taxonomic data, accession-to-taxon mappings, taxonomy branches, corrected ranks, and pathogen target taxa.

# Install

### Conda Formula
conda install bioconda::taxutils

### Pip Formula
pip install taxutils

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
tu = taxutils(accessions=None, low_memory=True, targets_json=None)

# Accession parsing and mapping
tu.parse_accession(header_strings, version=True)
tu.load_a2t(accessions, low_memory=None, extend=False)
tu.get_t2a(taxa, low_memory=None)

# Tree queries
tu.get_branch(taxon)
tu.get_subtree(taxon)
tu.get_lca(taxon_a, taxon_b)
tu.sort_taxa(taxa)

# Rank utilities
tu.get_rank_order()
tu.higher_than_rank(taxa, rank)
```

In taxutils, `accessions=list/of/accessions` can be passed to call load_a2t on construction of the taxutils object. A custom targets_json can similarly be passed in lieu of the default json explained below.  `load_a2t` overwrites `tu.a2t` by default; pass `extend=True` to add missing mappings without discarding existing ones. Method-level `low_memory=None` uses the mode set when `tu` was built.

# Rank correction

`taxutils` keeps the raw NCBI rank in `rank` and adds corrected rank columns. Canonical ranks (`R`, `D`, `K`, `P`, `C`, `O`, `F`, `G`, `S`) are used as anchors only when they move deeper than the corrected parent rank. Noncanonical ranks such as `no rank`, `clade`, and other unusual labels inherit position from the tree. If a child would be ranked at the same or a higher level than its parent, it is assigned a subrank such as `S2`, `S3`, or `F2`. The canonical name for the corrected rank is stored in `new_rank`.

# Target taxa

In ZarLab, we are working on metagenomics in the clinical setting, with the goal of creating an "agnostic diagnostic". We often want to look at broad array of taxa (`tu.target_taxa`) that could cause harm to people. In June 2024, CZI did the work of compiling a list of pathogenic taxa. I did the easy work of turning this into a json and uploading it to my website, so that it is available and easily accessed for all time (in case that link ever breaks). taxutils will extend the taxa list to include subtrees of each of those pathogenic taxa. It will additionally include SARS-CoV2, since it was excluded from CZI's list. If you find any other obvious, missing pathogens, please send me a note, so I can update my json. You can also update the target_taxa member variable yourself, or store an entirely different set of targets, if you wanted.

# Contact

Author: Will O'Brien  
Affiliation: Computer Science Department, UCLA  
Email: wob@cs.ucla.edu
