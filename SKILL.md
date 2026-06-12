---
name: taxutils
description: Use this skill when working with the Python taxutils package for NCBI taxonomy data, accession parsing from FASTA/Kraken-style headers, accession-to-taxid and taxid-to-accession lookups, taxonomic branches/subtrees/LCAs, corrected taxonomic ranks, pathogen target taxa, FASTA header cleaning, accession map generation, or Kraken-like taxonomy report formatting.
---

# taxutils

Use `taxutils` for local NCBI taxonomy workflows: parse accession IDs, map accessions to taxids, map taxids back to accessions, query the taxonomic tree, use corrected rank labels, and build pathogen-focused taxa sets.

## Import and Setup

Set `TAXUTILS_GLOBALS` before importing `taxutils` whenever the cache location matters. The package reads this environment variable at import time.

```python
import os
os.environ["TAXUTILS_GLOBALS"] = "/path/to/taxutils/cache"

from taxutils import taxutils
tu = taxutils()
```

If `TAXUTILS_GLOBALS` is not set, resources are stored under `./taxutils/` relative to the current working directory. Managed resources include `names.dmp`, `nodes.dmp`, `targets.json`, `nucl_gb.accession2taxid.gz`, and optionally `nucl_gb.accession2taxid.db`.

Use:

```python
tu = taxutils(accessions=None, low_memory=True, targets_json=None, rebuild=False)
```

- `accessions`: optional accession/header list to load into `tu.a2t` during construction.
- `low_memory=True`: default; scans compressed `nucl_gb.accession2taxid.gz` for requested lookups.
- `low_memory=False`: builds/reuses a SQLite accession database for faster repeated lookup work. Expect a slow first build and large disk usage.
- `targets_json`: custom pathogen/target JSON path in place of the default downloaded target list.
- `rebuild=True`: redownloads managed taxonomy/target/accession files and rebuilds the SQLite database.

## Core Object

The public constructor returns a `TaxonomicUtils` object:

```python
from taxutils import taxutils
tu = taxutils(low_memory=False)
```

Important members:

- `tu.names`: `{taxon: scientific_name}` including `0: "unclassified"` plus patched SARS-CoV-2 names.
- `tu.nodes`: pandas DataFrame with `taxon`, `parent`, raw NCBI `rank`, corrected `rank_code`, `rank_base`, `rank_idx`, and `new_rank`.
- `tu.parent`: `{taxon: parent_taxon}` with root `1` mapped to `None`.
- `tu.target_taxa`: taxonomically ordered pathogen target taxa.
- `tu.a2t`: accession-to-taxid map after construction with `accessions=` or after `load_a2t`.

The object repr lists available public methods.

## Accession Parsing

Use `tu.parse_accession` for raw accession strings, FASTA headers, Kraken-style headers, pandas Series, numpy arrays, or lists. Versions are kept by default.

```python
headers = [
    ">NC_045512.2 SARS-CoV-2",
    ">kraken:taxid|2886930|NC_001422.1 Escherichia phage phiX174",
]

acc_ids = tu.parse_accession(headers)                 # ["NC_045512.2", "NC_001422.1"]
acc_ids_no_version = tu.parse_accession(headers, version=False)
```

For FASTA files, read only header lines and pass them as a batch:

```python
headers = []
with open("input.fasta") as f:
    for line in f:
        if line.startswith(">"):
            headers.append(line.strip())

acc_ids = tu.parse_accession(headers)
```

When processing large FASTA files, stream headers in chunks instead of reading the whole file. Deduplicate accessions before lookup while preserving first-seen order:

```python
accessions = {}
for header in iter_header_lines("input.fasta"):
    for accession in tu.parse_accession(header):
        accessions.setdefault(accession, None)
acc_ids = list(accessions)
```

## Accession and Taxid Maps

Load an accession-to-taxid subset into `tu.a2t`:

```python
tu.load_a2t(acc_ids)
taxon = tu.a2t[acc_ids[0]]
name = tu.names[taxon]
```

`load_a2t` parses input strings with versions enabled and overwrites `tu.a2t` by default. Preserve existing mappings with:

```python
tu.load_a2t(more_acc_ids, extend=True)
```

Use `get_t2a` for direct taxid-to-accession lookup. It returns accessions assigned to the provided taxa only; pass a subtree to include descendants.

```python
enterovirus = 12059
direct_accessions = tu.get_t2a([enterovirus])
subtree_accessions = tu.get_t2a(tu.get_subtree(enterovirus))
```

Prefer `low_memory=False` for repeated `load_a2t` or `get_t2a` calls in notebooks, pipelines, or scripts that can afford the SQLite database. Use default low-memory mode for one-off lookups or constrained disk environments.

## FASTA Workflows

To create a two-column accession-to-taxid map for external tools:

```python
tu = taxutils(low_memory=False)

seen = {}
with open("input.fasta") as in_f:
    for line in in_f:
        if line.startswith(">"):
            for accession in tu.parse_accession(line):
                seen.setdefault(accession, None)

acc_ids = list(seen)
tu.load_a2t(acc_ids)

with open("taxid_map.tsv", "w") as f:
    for accession in acc_ids:
        taxon = tu.a2t.get(accession)
        if taxon is not None:
            f.write(f"{accession}\t{taxon}\n")
```

To clean FASTA headers to accession-only headers:

```python
accessions = tu.parse_accession(header_line)
if not accessions:
    raise ValueError(f"No accession found: {header_line.strip()}")
clean_header = f">{accessions[0]}\n"
```

For in-place FASTA rewrites, write to a temporary file in the same directory and then use `os.replace(tmp_path, input_path)` so partial writes do not corrupt the input.

## Tree Queries

Use these public methods:

```python
branch = tu.get_branch(taxon)      # root-to-taxon branch
subtree = tu.get_subtree(taxon)    # taxon plus descendants
lca = tu.get_lca(taxon_a, taxon_b)
distance = tu.get_distance(taxon_a, taxon_b)
ordered = tu.sort_taxa(taxa)
tree = tu.format_tree(taxa)        # Series indexed by taxon with indented names
```

`format_tree(taxa, include_ancestors=True, root=1, indent="\t")` includes ancestors by default and returns a pandas Series named `name`.

For branch rows in branch order:

```python
branch = tu.get_branch(taxon)
branch_df = (
    tu.nodes.loc[tu.nodes["taxon"].isin(branch)]
    .set_index("taxon")
    .loc[branch]
    .reset_index()
)
```

For comparing two taxa:

```python
branch_a = tu.get_branch(taxon_a)
branch_b = tu.get_branch(taxon_b)
lca = tu.get_lca(taxon_a, taxon_b)
distance = tu.get_distance(taxon_a, taxon_b)
branch_taxa = tu.sort_taxa(set(branch_a) | set(branch_b))
```

Use `tu.nodes` and `tu.names` to resolve unknown taxa by name instead of relying on memorized taxids:

```python
tu.nodes["name"] = tu.nodes["taxon"].map(tu.names)
hits = tu.nodes[tu.nodes["name"].str.lower().str.contains("rhinovirus c22", na=False)]
```

## Rank Utilities

Raw NCBI rank is preserved in `rank`. Corrected rank columns are:

- `rank_code`: corrected rank code, including subrank suffixes such as `F2`, `G2`, or `S3`.
- `rank_base`: base canonical code.
- `rank_idx`: numeric rank order.
- `new_rank`: canonical corrected rank name.

Canonical order:

```python
tu.get_rank_order()
# ["U", "R", "D", "K", "P", "C", "O", "F", "G", "S"]
```

Rank aliases accepted by `higher_than_rank` include codes and names such as `F`, `family`, `subfamily`, `D`, `domain`, `superkingdom`, `C`, and `clade`.

```python
mask = tu.higher_than_rank(tu.nodes["taxon"], "F")
at_or_below_family = ~mask
```

`higher_than_rank` returns `True` for taxa higher than the threshold and `False` for taxa at or below it. For unknown taxids, it defaults to the threshold and therefore returns `False`.

To inspect NCBI/corrected rank mismatches:

```python
tu.nodes["name"] = tu.nodes["taxon"].map(tu.names)
mismatches = tu.nodes[tu.nodes["new_rank"] != tu.nodes["rank"]]
```

## Target Taxa

`tu.target_taxa` is built from `targets_json`, defaulting to the package-managed pathogen target list. The package expands target taxa through descendants and relevant lower-rank ancestors, then returns them in taxonomic order.

Use it directly for pathogen filtering:

```python
target_species = tu.nodes[
    tu.nodes["taxon"].isin(tu.target_taxa)
    & tu.nodes["rank_base"].eq("S")
]
```

For a custom target universe, pass `targets_json=` during construction or assign your own ordered list/set to `tu.target_taxa`.

## Kraken Read-Level Movement Analysis

For Kraken read-level classification output, use `taxutils(low_memory=False)`, parse the accession column in one batch, and use the package target set directly:

```python
from taxutils import taxutils

tu = taxutils(low_memory=False)
out = pd.read_csv(
    "kraken_output.tsv",
    sep="\t",
    header=None,
    names=["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"],
)

out["accession"] = tu.parse_accession(out["accession"])
tu.load_a2t(out["accession"].dropna().unique().tolist())
out["labeled_taxid"] = out["accession"].map(tu.a2t)

targets = set(tu.target_taxa)
out["is_target"] = out["labeled_taxid"].isin(targets)
```

For each taxon observed in a Kraken `lca_mapping`, compare it to the known labeled taxon with direct tree calls:

```python
lca_taxon = tu.get_lca(labeled_taxon, observed_taxon)
distance = tu.get_distance(labeled_taxon, observed_taxon)

if observed_taxon in targets and lca_taxon not in targets:
    movement = f"{observed_taxon}->{lca_taxon}"
```

For repeated rank threshold checks in a row loop, build lookup dictionaries from `tu.nodes` once:

```python
rank_idx = tu.nodes.set_index("taxon")["rank_idx"].to_dict()
rank_order = {rank: i for i, rank in enumerate(tu.get_rank_order())}

if rank_idx.get(lca_taxon, rank_order["F"]) < rank_order["F"]:
    common_family_level_count += count
```

## Kraken-like Report Formatting

Use `format_tree` plus `get_branch` to build count reports similar to Kraken2 `kreport` output.

```python
taxon_counts = {12059: 500, 147711: 120, 2697049: 80}
tree = tu.format_tree(taxon_counts.keys())

cumulative = {taxon: 0 for taxon in tree.index}
for taxon, count in taxon_counts.items():
    for parent_taxon in tu.get_branch(taxon):
        if parent_taxon in cumulative:
            cumulative[parent_taxon] += count

total = sum(taxon_counts.values())
ranks = tu.nodes.set_index("taxon")[["rank_code", "new_rank", "rank"]]

df = tree.to_frame("name")
df["pct"] = [round(100 * cumulative[t] / total, 2) if total else 0 for t in df.index]
df["cumulative_count"] = [cumulative[t] for t in df.index]
df["count"] = [taxon_counts.get(t, 0) for t in df.index]
df = df.join(ranks).fillna({"rank_code": "U", "new_rank": "U", "rank": "U"})
df["taxon"] = df.index
df = df[["pct", "cumulative_count", "count", "rank_code", "new_rank", "rank", "taxon", "name"]]
```

Use `tu.get_lca(a, b)` to verify that reported taxa preserve expected hierarchy.

## Practical Guidance

- Use `load_a2t` for accession subsets and `get_t2a` for selected taxon-to-accession lookups.
- Keep accession versions when mapping against NCBI `nucl_gb.accession2taxid.gz`; the package’s lookup key is `accession.version`.
- Use `set.update(...)` when adding many branch or subtree taxa to a set.
- Call `tu.sort_taxa(...)` after set operations whenever display order matters.
- Add `tu.nodes["name"] = tu.nodes["taxon"].map(tu.names)` before vectorized name searches or report tables.
- Use `rebuild=True` only when intentionally refreshing downloaded taxonomy/accession/target resources.
