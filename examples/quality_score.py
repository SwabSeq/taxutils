import argparse
import os
from collections import OrderedDict

import pandas as pd

from taxutils import taxutils

KRAKEN_COLUMNS = ["status", "accession", "predicted_taxid", "seqlen", "lca_mapping"]
LABEL_COLUMNS = ["status", "accession", "label_taxon", "seqlen", "lca_mapping", "a2t_taxon"]
OUTPUT_COLUMNS = [
    "accession",
    "label_taxon",
    "a2t_taxon",
    "seq_len",
    "total_kmers",
    "unclassified_kmers",
    "classified_kmers",
    "masked_kmers",
    "total_distance_moved",
    "average_distance_moved",
    "topology_scale",
    "total_normalized_distance",
    "average_normalized_distance",
    "total_upward_loss",
    "average_upward_loss",
    "total_lateral_shift",
    "average_lateral_shift",
    "distance_variance",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate taxonomy distance quality metrics from Kraken output."
    )
    parser.add_argument(
        "-i",
        "--kraken-output",
        required=True,
        help="Path to Kraken standard read-level classification output.",
    )
    parser.add_argument(
        "-o",
        "--results",
        required=True,
        help="Path to output CSV.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print Kraken accession strings where no accession is found.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        default=50_000,
        help="Number of input rows to process per chunk.",
    )
    parser.add_argument(
        "--lca-cache-size",
        type=int,
        default=20_000,
        help="Deprecated; accepted for compatibility but has no effect. Raw mappings are no longer cached.",
    )
    parser.add_argument(
        "--lca-fill-cache-size",
        type=int,
        default=100_000,
        help="Maximum (accession, a2t_taxon) mappings to retain when filling blank lca_mapping rows.",
    )
    parser.add_argument(
        "--movement-cache-size",
        type=int,
        default=250_000,
        help="Maximum (label_taxon, observed_taxon) movement calculations to cache. Use 0 to disable.",
    )
    return parser.parse_args()


def aggregate_lca_mapping(lca_mapping):
    counts = {}
    for item in str(lca_mapping).split():
        # Kraken emits mate-pair and translated reading-frame boundaries.
        if item in {"|:|", "-:-"} or ":" not in item:
            continue
        taxon, count = item.split(":", 1)
        taxon = -1 if taxon == "A" else int(taxon)
        counts[taxon] = counts.get(taxon, 0) + int(count)
    return counts


def _remember(cache, key, value, max_size):
    if max_size <= 0:
        return
    cache[key] = value
    cache.move_to_end(key)
    while len(cache) > max_size:
        cache.popitem(last=False)


def infer_input_columns(path):
    first_row = pd.read_csv(path, sep="\t", header=None, nrows=1)
    n_columns = first_row.shape[1]
    if n_columns == 5:
        return KRAKEN_COLUMNS
    if n_columns == 6:
        return LABEL_COLUMNS
    raise ValueError(
        "Input must have 5 Kraken columns or 6 label-enumeration columns."
    )


def read_label_chunks(path, tu, chunksize, verbose=False):
    columns = infer_input_columns(path)
    for out in pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=columns,
        chunksize=chunksize,
    ):
        parsed_accessions = tu.parse_accession(out["accession"])
        if verbose:
            for value in out.loc[parsed_accessions == "NA", "accession"]:
                print(f"NA accession: {value}")
        out["accession"] = parsed_accessions
        out = out[out["accession"] != "NA"].copy()
        if out.empty:
            continue

        if columns == KRAKEN_COLUMNS:
            tu.load_a2t(out["accession"].drop_duplicates().tolist(), extend=True)
            out["label_taxon"] = pd.to_numeric(
                out["accession"].map(tu.a2t),
                errors="coerce",
            )
            out["a2t_taxon"] = out["label_taxon"]
            yield out.loc[:, LABEL_COLUMNS].copy()
            continue

        out["lca_mapping"] = out["lca_mapping"].replace({
            "": pd.NA,
            "None": pd.NA,
            "nan": pd.NA,
        })
        yield out.loc[:, LABEL_COLUMNS].copy()


def fill_lca_mappings(out, lca_fill_state, max_size):
    """Parse mappings and fill missing rows using cached taxon counts."""
    raw_mappings = out["lca_mapping"]
    out = out.drop(columns="lca_mapping")
    cache = lca_fill_state["cache"]
    filled = []
    for accession, a2t_taxon, lca_mapping in zip(
        out["accession"],
        out["a2t_taxon"],
        raw_mappings,
    ):
        key = None if pd.isna(a2t_taxon) else (accession, int(a2t_taxon))
        if pd.notna(lca_mapping):
            lca_mapping = aggregate_lca_mapping(lca_mapping)
            lca_fill_state["last_key"] = key
            lca_fill_state["last_lca_counts"] = lca_mapping
            if key is not None and max_size > 0:
                _remember(cache, key, lca_mapping, max_size)
            filled.append(lca_mapping)
        elif key is not None and key == lca_fill_state["last_key"]:
            filled.append(lca_fill_state["last_lca_counts"])
        elif max_size > 0 and key in cache:
            cache.move_to_end(key)
            filled.append(cache[key])
        else:
            filled.append(pd.NA)
    # Keep variable-sized mappings as Python objects. Assigning a list of
    # strings directly makes NumPy infer a fixed-width Unicode dtype and can
    # attempt an enormous temporary allocation for long Kraken mappings.
    out["lca_counts"] = pd.Series(filled, index=out.index, dtype=object)
    return out


def normalize_chunk(out, lca_fill_state, lca_fill_cache_size):
    out = out.copy()
    taxon_columns = ["label_taxon", "a2t_taxon"]
    numeric_taxa = out.loc[:, taxon_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    out.loc[:, taxon_columns] = numeric_taxa
    out = fill_lca_mappings(out, lca_fill_state, lca_fill_cache_size)
    out = out.dropna(subset=["label_taxon", "a2t_taxon", "lca_counts"]).copy()
    if out.empty:
        return out
    return out.astype({"label_taxon": int, "a2t_taxon": int})


def update_topology_scales(tu, label_taxa, topology_scales):
    missing = [
        int(taxon)
        for taxon in pd.unique(label_taxa)
        if int(taxon) not in topology_scales
    ]
    if not missing:
        return
    topology_scales.update(
        tu.topology(
            missing,
            anchor_rank="F",
            stat="topology_scale",
        ).to_dict()
    )


def movement_stats(tu, labeled_taxon, taxon, cache, max_size):
    key = (int(labeled_taxon), int(taxon))
    if max_size > 0 and key in cache:
        cache.move_to_end(key)
        return cache[key]

    lca_taxon = tu.get_lca(labeled_taxon, taxon)
    lca_depth = tu._get_depth(lca_taxon)
    upward_loss = tu._get_depth(labeled_taxon) - lca_depth
    lateral_shift = tu._get_depth(taxon) - lca_depth
    distance = upward_loss + lateral_shift
    stats = (upward_loss, lateral_shift, distance)
    _remember(cache, key, stats, max_size)
    return stats


def write_results_chunk(path, rows, include_header):
    if not rows:
        return include_header

    results = pd.DataFrame.from_records(rows, columns=OUTPUT_COLUMNS)
    results.to_csv(
        path,
        mode="w" if include_header else "a",
        header=include_header,
        index=False,
    )
    return False


def write_empty_results(path):
    pd.DataFrame(columns=OUTPUT_COLUMNS).to_csv(path, index=False)


def distance_metrics(
    tu,
    labeled_taxon,
    kmer_counts,
    topology_scale,
    movement_cache=None,
    movement_cache_size=0,
):
    kmer_counts = dict(kmer_counts)
    total_kmers = sum(kmer_counts.values())
    unclassified_kmers = kmer_counts.pop(0, 0)
    masked_kmers = kmer_counts.pop(-1, 0)
    classified_kmer_count = sum(kmer_counts.values())

    total_distance_moved = 0
    total_upward_loss = 0
    total_lateral_shift = 0
    distance_sq_sum = 0

    for taxon, count in kmer_counts.items():
        upward_loss, lateral_shift, distance = movement_stats(
            tu,
            labeled_taxon,
            taxon,
            movement_cache,
            movement_cache_size,
        )

        total_distance_moved += count * distance
        total_upward_loss += count * upward_loss
        total_lateral_shift += count * lateral_shift
        distance_sq_sum += count * distance * distance

    average_distance_moved = (
        total_distance_moved / classified_kmer_count
        if classified_kmer_count
        else 0
    )
    average_upward_loss = (
        total_upward_loss / classified_kmer_count if classified_kmer_count else 0
    )
    average_lateral_shift = (
        total_lateral_shift / classified_kmer_count if classified_kmer_count else 0
    )
    total_normalized_distance = total_distance_moved / topology_scale
    average_normalized_distance = (
        total_normalized_distance / classified_kmer_count
        if classified_kmer_count
        else 0
    )
    distance_variance = (
        distance_sq_sum / classified_kmer_count - average_distance_moved ** 2
        if classified_kmer_count
        else 0
    )

    return {
        "total_kmers": total_kmers,
        "unclassified_kmers": unclassified_kmers,
        "classified_kmers": classified_kmer_count,
        "masked_kmers": masked_kmers,
        "total_distance_moved": total_distance_moved,
        "average_distance_moved": average_distance_moved,
        "topology_scale": topology_scale,
        "total_normalized_distance": total_normalized_distance,
        "average_normalized_distance": average_normalized_distance,
        "total_upward_loss": total_upward_loss,
        "average_upward_loss": average_upward_loss,
        "total_lateral_shift": total_lateral_shift,
        "average_lateral_shift": average_lateral_shift,
        "distance_variance": max(distance_variance, 0),
    }


def main():
    args = parse_args()
    tu = taxutils(low_memory=False, keep_accession_downloads=False)

    results_dir = os.path.dirname(args.results)
    if results_dir:
        os.makedirs(results_dir, exist_ok=True)

    topology_scales = {}
    movement_cache = OrderedDict()
    lca_fill_state = {
        "cache": OrderedDict(),
        "last_key": None,
        "last_lca_counts": pd.NA,
    }
    include_header = True

    for chunk in read_label_chunks(
        args.kraken_output,
        tu,
        chunksize=args.chunksize,
        verbose=args.verbose,
    ):
        chunk = normalize_chunk(
            chunk,
            lca_fill_state=lca_fill_state,
            lca_fill_cache_size=args.lca_fill_cache_size,
        )
        if chunk.empty:
            continue

        update_topology_scales(tu, chunk["label_taxon"], topology_scales)

        rows = []
        score_columns = [
            "accession",
            "label_taxon",
            "a2t_taxon",
            "seqlen",
            "lca_counts",
        ]
        for accession, label_taxon, a2t_taxon, seqlen, lca_counts in chunk[
            score_columns
        ].itertuples(index=False, name=None):
            metrics = distance_metrics(
                tu=tu,
                labeled_taxon=label_taxon,
                kmer_counts=lca_counts,
                topology_scale=topology_scales[label_taxon],
                movement_cache=movement_cache,
                movement_cache_size=args.movement_cache_size,
            )
            rows.append(
                {
                    "accession": accession,
                    "label_taxon": label_taxon,
                    "a2t_taxon": a2t_taxon,
                    "seq_len": seqlen,
                    **metrics,
                }
            )

        include_header = write_results_chunk(args.results, rows, include_header)

    if include_header:
        write_empty_results(args.results)

if __name__ == "__main__":
    main()
