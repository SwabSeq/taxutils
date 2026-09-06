---
name: taxutils
description: Use this skill when working with the Python taxutils package for NCBI taxonomy data, accession parsing from FASTA/Kraken-style headers, accession-to-taxid and taxid-to-accession lookups, taxonomic branches/subtrees/LCAs, corrected taxonomic ranks, pathogen target taxa, FASTA header cleaning, accession map generation, or Kraken-like taxonomy report formatting.
---

# taxutils

Use `taxutils` for local NCBI taxonomy workflows: parse accession IDs, map accessions to taxids, map taxids back to accessions, query the taxonomic tree, use corrected rank labels, and build pathogen-focused taxa sets.

## Import and Setup

Pass `save_folder=` whenever the cache location matters. `TAXUTILS_GLOBALS` is the default when `save_folder` is not given, and is the only environment variable the package reads; it is read on each call, so setting it after import works.

```python
from taxutils import backend_info, taxutils

tu = taxutils(save_folder="/path/to/taxutils/cache")
print(backend_info())
```

`threads` sets the worker count for every parallel stage, including the accession database build; `None` uses all logical CPUs.

```python
tu = taxutils(low_memory=False, threads=8)
```

If `TAXUTILS_GLOBALS` is not set, resources are stored under `./taxutils/` relative to the current working directory. Managed resources include `names.dmp`, `nodes.dmp`, `targets.json`, `nucl_gb.accession2taxid.gz`, optionally `nucl_wgs.accession2taxid.gz`, and optionally `nucl.accession2taxid.db`.

Use:

```python
tu = taxutils(
    accessions=None,
    low_memory=True,
    targets_json=None,
    wgs=False,
    save_folder=None,
    threads=None,
    refresh=False,
)
```

- `accessions`: optional accession/header list to load into `tu.a2t` during construction.
- `low_memory=True`: default; scans compressed accession2taxid files for requested lookups.
- `low_memory=False`: builds/reuses a SQLite accession database for faster repeated lookup work. The first build is the expensive one; it decompresses the sources in parallel and merges them into the table in key order.
- `targets_json`: custom pathogen/target JSON path in place of the default downloaded target list.
- `refresh=True`: re-fetches the managed taxonomy/target files and brings an existing SQLite database up to date, applying only the accessions NCBI added, changed, or withdrew. Skips the work entirely when the sources are unchanged. There is no separate rebuild switch: anything missing is downloaded, and a missing or unreadable database is rebuilt automatically.
- `wgs=False`: default; uses `nucl_gb.accession2taxid.gz` only. Pass `wgs=True` to also download/use `nucl_wgs.accession2taxid.gz` for WGS/TSA accessions.

SQLite mode always uses `nucl.accession2taxid.db`; if it was built GB-only, a later `wgs=True` call upgrades the same DB with WGS mappings. A database written before the current schema is detected and rebuilt once, automatically.

Long backend calls are interruptible: `Ctrl-C` during a build or lookup raises `KeyboardInterrupt` promptly, discards the partial database, and leaves any installed database untouched.

## Backend and Integration Strategy

The native Rust extension is mandatory. Published Python wheels bundle it, and
imports fail clearly if it is missing or exposes an incompatible API. There is
no `TAXUTILS_BACKEND` mode switch and no pure-Python fallback. Use
`backend_info()` for diagnostics:

```python
from taxutils import backend_info
assert backend_info()["selected"] == "rust"
```

Source builds require a Rust toolchain but not a sibling repository checkout.
The extension depends on `taxutils >=1.1.1,<2` from crates.io. Refresh
`rust/Cargo.lock` after compatible crate releases so Python wheels inherit the
new backend; do not copy the crate source into this repository.

Call the public Python APIs. Do not import `taxutils._rust`, shell out to the
`tu` executable from Python, or copy FASTA, accession lookup, download, or
database-building implementations into custom scripts. The public wrappers
provide compatibility checks, stable errors and return types, and the fastest
maintained implementation without creating redundant code.

Rust currently owns SQLite database construction, bulk accession lookups, and
the extract/clean/grep/filter FASTA engines exclusively. Taxonomy methods that
naturally operate on pandas or NumPy containers stay in Python so they do not
pay unnecessary conversion costs.

When adding a custom Python workflow, compose the public functions in this
package and batch work across the Python/Rust boundary. Add a new public Rust
crate API plus a thin PyO3/Python wrapper only when existing operations cannot
express the workflow. Never add a second Python implementation as a fallback.
When adding a custom Rust binary or script, depend on the `taxutils` crate and
call its library APIs rather than invoking the `tu` CLI as a subprocess.

## Choosing Lookup Mode and Disk Layout

- Use `low_memory=True` for a one-off lookup, constrained storage, or workflows
  that already retain the compressed NCBI mappings. Each lookup scans gzip
  input.
- Use `low_memory=False` for repeated lookups, reverse taxid-to-accession
  queries, FASTA filtering, notebooks, and pipelines. The Rust backend builds
  the shared SQLite index atomically and reuses it.
- The compressed NCBI inputs are always retained, because low-memory lookups
  scan them directly. On a SQLite-only cache they can be deleted by hand to
  reclaim the space; a later low-memory lookup downloads them again.
- Set `wgs=True` only when WGS/TSA accessions are needed; it materially expands
  download, build, and disk requirements.

Build the indexed database once in a persistent save folder and reuse that
cache across scripts. Use `refresh=True` only when new NCBI data is wanted; a
plain call reuses whatever is already installed without touching the network.

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

Use `tu.parse_accession` for raw accession strings, FASTA headers, Kraken-style headers, pandas Series, numpy arrays, or lists. Versions are kept by default. It returns the first accession found from each string using the same container type where possible; missing accessions are returned as `"NA"`.

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
    accession = tu.parse_accession(header)
    if accession != "NA":
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

`load_a2t` parses input strings with versions enabled and overwrites `tu.a2t` by default. Method-level `wgs=None` uses the constructor setting; pass `wgs=True` to include WGS/TSA accessions for a specific call. Preserve existing mappings with:

```python
tu.load_a2t(more_acc_ids, extend=True)
```

Use `get_t2a` for direct taxid-to-accession lookup. It returns accessions assigned to the provided taxa only; pass a subtree to include descendants.

```python
enterovirus = 12059
direct_accessions = tu.get_t2a([enterovirus])
subtree_accessions = tu.get_t2a(tu.get_subtree(enterovirus))
```

Prefer `low_memory=False` for repeated `load_a2t` or `get_t2a` calls in notebooks, pipelines, or scripts that can afford the SQLite database. Use default low-memory mode for one-off lookups or constrained disk environments. Pass deduplicated accessions in batches rather than calling `load_a2t` once per row.

## FASTA Workflows

Call the package functions directly; each is a thin wrapper around the required
Rust engine:

```python
from taxutils.clean_headers import clean_fasta_headers
from taxutils.extract_headers import extract_accessions
from taxutils.filter_fasta import filter_fasta
from taxutils.grep_fasta import grep_fasta

extract_accessions("input.fasta", "accessions.txt")
clean_fasta_headers("input.fasta", "clean.fasta")
grep_fasta("input.fasta", "NC_045512.2", "hits.fasta")
filter_fasta("input.fasta", "viral.fasta", {2697049}, filter_mode="keep")
```

These functions already stream in bounded batches, preserve record order, and
use atomic output where applicable. Do not recreate their record loops in
examples or custom scripts unless the workflow needs genuinely different
semantics.

To create a two-column accession-to-taxid map for external tools:

```python
tu = taxutils(low_memory=False)

seen = {}
with open("input.fasta") as in_f:
    for line in in_f:
        if line.startswith(">"):
            accession = tu.parse_accession(line)
            if accession != "NA":
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
accession = tu.parse_accession(header_line)
if accession == "NA":
    raise ValueError(f"No accession found: {header_line.strip()}")
clean_header = f">{accession}\n"
```

For in-place FASTA rewrites, write to a temporary file in the same directory and then use `os.replace(tmp_path, input_path)` so partial writes do not corrupt the input.

## Tree Queries

Use these public methods:

```python
branch = tu.get_branch(taxon)      # root-to-taxon branch
subtree = tu.get_subtree(taxon)    # taxon plus descendants
family = tu.get_ancestor(taxon, anchor_rank="F")
leaf = tu.is_leaf(taxon)           # True if taxon has no child nodes
child = tu.is_child(taxon_a, taxon_b)
descendent = tu.is_descendent(taxon_a, taxon_b)
lca = tu.get_lca(taxon_a, taxon_b)
distance = tu.get_distance(taxon_a, taxon_b)
ordered = tu.sort_taxa(taxa)
tree = tu.format_tree(taxa)        # Series indexed by taxon with indented names
profile = tu.topology(taxon)       # Series of subtree topology metrics
scale = tu.topology(taxon, anchor_rank="F", stat="topology_scale")
```

`format_tree(taxa, include_ancestors=True, root=1, indent="\t")` includes ancestors by default and returns a pandas Series named `name`.

`get_ancestor(taxon, anchor_rank)` returns the nearest ancestor at the requested corrected rank. If the input taxon already has that rank, it returns the input taxon; if no ancestor has that rank, it returns the input taxon as a fallback. It accepts a scalar taxon, list-like input, numpy arrays, or pandas Series and returns the same container type where possible.

`is_leaf(taxon)` returns whether taxa have no children in the taxonomy tree. A scalar taxon returns a `bool`; a list-like input returns a list of booleans; a numpy array returns a boolean array with the original shape; and a pandas Series returns a boolean Series with the original index.

`is_child(taxon_a, taxon_b)` returns whether `taxon_a` is a direct child of `taxon_b`. It uses `tu.parent` directly, so each pairwise check is O(1).

`is_descendent(taxon_a, taxon_b)` returns whether `taxon_a` is a strict descendant of `taxon_b`; a taxon is not considered a descendant of itself. The first call builds a cached tree interval index in O(n), then each pairwise check is O(1).

`is_child` and `is_descendent` accept either two scalar taxa or two list-like inputs of the same length. Scalar inputs return a `bool`; list-like inputs return a list of booleans; numpy arrays return boolean arrays with the original `taxon_a` shape; and pandas Series return boolean Series with the original `taxon_a` index.

`topology(taxon, anchor_rank=None, stat=None)` returns subtree topology metrics including taxon counts, leaf fraction, depth, branchiness, and `topology_scale`. Pass `anchor_rank="F"` to summarize the nearest family ancestor instead of the exact taxon. With `stat=None`, a single taxon returns a pandas Series and a list, array, or pandas Series returns a DataFrame.

Pass `stat` to return one metric. Valid values are `n_taxa`, `n_leaves`, `max_depth`, `mean_depth`, `topology_scale`, `max_children`, `branching_taxa_fraction`, and `top_child_fraction`. A scalar taxon returns a scalar; a list, array, or pandas Series returns a Series indexed by taxon.

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

For Kraken read-level classification output, use indexed mode, parse the accession column directly, and use the package target set:

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

## Native Rust Scripts

When the custom program itself is Rust, depend on the `taxutils` crate and use
its library API rather than invoking the CLI. Keeping a compatible semver range
allows package updates to pick up backend improvements after the lockfile is
refreshed.

```toml
[dependencies]
taxutils = ">=1.1.1, <2"
```

For explicit database preparation:

```rust
use taxutils::{ensure_accession_database, AccessionDatabaseOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let database = ensure_accession_database(
        "/path/to/taxutils/cache",
        AccessionDatabaseOptions {
            wgs: false,
            keep_downloads: false,
            ..Default::default()
        },
    )?;
    println!("{}", database.display());
    Ok(())
}
```

Use the crate's `lookup_accession_taxids`, `lookup_taxid_accessions`,
`extract_accessions`, `clean_fasta_headers`, `grep_fasta`, and `filter_fasta`
APIs for custom Rust workflows. Database construction and upgrades are built in
a sibling temporary file and installed only after all rows and indexes are
complete.

## Practical Guidance

- Use `load_a2t` for accession subsets and `get_t2a` for selected taxon-to-accession lookups.
- Check `backend_info()` when performance is operationally important.
- Import public FASTA functions from their `taxutils.*` modules; never call the private native extension directly.
- Deduplicate and batch accession lookups; avoid per-row database calls.
- Keep accession versions when mapping against NCBI accession2taxid files; the package’s lookup key is `accession.version`.
- Use `set.update(...)` when adding many branch or subtree taxa to a set.
- Call `tu.sort_taxa(...)` after set operations whenever display order matters.
- Add `tu.nodes["name"] = tu.nodes["taxon"].map(tu.names)` before vectorized name searches or report tables.
- Use `refresh=True` to pick up new NCBI data; it updates the installed database in place instead of discarding it.
