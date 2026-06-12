---
name: taxutils
description: Use this skill when working with the Python taxutils package for NCBI taxonomy data, accession parsing, accession-to-taxon and taxon-to-accession lookups, taxonomy branches/subtrees/LCA, corrected ranks, or pathogen target taxa.
---

# taxutils Usage Skill

## Public entry point

Use the package constructor from the top-level package:

```python
from taxutils import taxutils

tu = taxutils(low_memory=False)
```

## Setup

Use `TAXUTILS_GLOBALS` to control where taxutils stores downloaded taxonomy files, target metadata, accession maps, and the optional SQLite database. Recommend that users export one consistent `TAXUTILS_GLOBALS` folder in their shell config. For a temporary script or notebook, set `os.environ` before importing taxutils:

```python
import os
os.environ["TAXUTILS_GLOBALS"] = "/path/to/taxutils/saves"

from taxutils import taxutils
tu = taxutils()
```

Default mode is low-memory and scans the compressed accession map when needed:

```python
tu = taxutils()
```

For faster repeated accession lookups, build/reuse the SQLite database:

```python
tu = taxutils(low_memory=False)
```

To refresh managed resource files and rebuild the accession SQLite database:

```python
tu = taxutils(rebuild=True)
```

## Core object

Important public members:

- `tu.names`: `{taxon: scientific_name}`
- `tu.nodes`: taxonomy dataframe with `taxon`, `parent`, raw `rank`, corrected `rank_code`, `rank_base`, `rank_idx`, and `new_rank`
- `tu.target_taxa`: pathogen-derived target taxa
- `tu.a2t`: accession-to-taxon map after `load_a2t`
- `tu.parent`: `{taxon: parent_taxon}`

The repr lists public methods.

## Accession parsing

Use the instance method:

```python
accessions = tu.parse_accession(headers_or_strings)
accessions_no_version = tu.parse_accession(headers_or_strings, version=False)
```

Inputs may be a string, list, pandas Series, or numpy array. Strings can be raw accession IDs or larger FASTA headers. Keep versions by default.

For a FASTA file:

```python
headers = []
with open("input.fasta") as f:
    for line in f:
        if line.startswith(">"):
            headers.append(line.strip())

acc_ids = tu.parse_accession(headers)
```

## Accession-to-taxon and taxon-to-accession

Load a subset accession map:

```python
tu.load_a2t(acc_ids)
taxon = tu.a2t[acc_ids[0]]
name = tu.names[taxon]
```

`load_a2t` overwrites `tu.a2t` by default. To preserve the existing map and load only missing accessions:

```python
tu.load_a2t(more_acc_ids, extend=True)
```

Get accessions for taxa with `get_t2a`:

```python
accessions = tu.get_t2a([12059])
subtree_accessions = tu.get_t2a(tu.get_subtree(12059))
```

## Tree queries

Use:

```python
branch = tu.get_branch(taxon)      # root-to-taxon branch
subtree = tu.get_subtree(taxon)    # taxon plus descendants
lca = tu.get_lca(taxon_a, taxon_b)
ordered = tu.sort_taxa(taxa)
tree = tu.format_tree(taxa)        # Series indexed by taxon with indented names
```

When displaying branch rows in branch order:

```python
branch = tu.get_branch(taxon)
tu.nodes.loc[tu.nodes["taxon"].isin(branch)].set_index("taxon").loc[branch].reset_index()
```

To combine branches:

```python
branch_a = tu.get_branch(taxon_a)
branch_b = tu.get_branch(taxon_b)
branch_taxa = tu.sort_taxa(set(branch_a) | set(branch_b))
```

To format taxa as an indented tree:

```python
import pandas as pd

taxon_counts = pd.Series({12059: 500, 147711: 120}, name="count")
tree = tu.format_tree(taxon_counts.index)

report = pd.DataFrame({
    "count": taxon_counts.reindex(tree.index).fillna(0).astype(int),
    "name": tree,
})

with open("example_tree.txt", "w") as f:
    for taxon, row in report.iterrows():
        f.write(f"{taxon}\t{row['count']}\t{row['name']}\n")
```

## Rank correction

Raw NCBI rank is kept in `rank`. Corrected columns are:

- `rank_code`: e.g. `F`, `F2`, `S`, `S2`
- `rank_base`: canonical base code
- `rank_idx`: canonical rank order index
- `new_rank`: canonical rank name

Canonical rank order:

```python
tu.get_rank_order()
# ['U', 'R', 'D', 'K', 'P', 'C', 'O', 'F', 'G', 'S'] 
```

The rank order codes can be mapped in nodes from `rank_base` to new_rank, rank, etc. "U" is unclassified.


Rank filtering:

```python
mask = tu.higher_than_rank(tu.nodes["taxon"], "F")
```

`higher_than_rank` returns `True` for taxa higher than the passed rank and `False` for taxa at or below it.

## Pandas patterns

Add names to nodes when needed:

```python
tu.nodes["name"] = tu.nodes["taxon"].map(tu.names)
```

Search names safely:

```python
hits = tu.nodes[tu.nodes["name"].str.lower().str.contains("rhinovirus c22", na=False)]
```

Use `set.update(...)` to add many taxa to a set:

```python
taxa = set(tu.get_branch(taxon))
taxa.update(tu.get_subtree(other_taxon))
```
