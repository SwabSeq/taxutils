"""Offline regression tests: python -m unittest discover -s tests -p 'test_quality_score.py'."""

import argparse
from collections import OrderedDict
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "examples"))

import quality_score as qs
from taxutils.taxutils import TaxonomicUtils


def fill_state():
    return {"cache": OrderedDict(), "last_key": None, "last_lca_counts": pd.NA}


def fixture_taxonomy():
    # Family 10 contains genera 11/12, species 13/14 and 15, and strain 16.
    parents = {1: None, 10: 1, 11: 10, 12: 10, 13: 11, 14: 11, 15: 12, 16: 13}
    nodes = pd.DataFrame({
        "taxon": list(parents),
        "parent": list(parents.values()),
        "rank_base": ["R", "F", "G", "G", "S", "S", "S", "S"],
    })
    return TaxonomicUtils({}, nodes, [], a2t={"NC_000001.1": 13}, parent=parents)


class QualityScoreTests(unittest.TestCase):
    mapping = "13:2 0:5 11:2 15:4 A:2 |:| 14:1 -:- 13:1"

    def setUp(self):
        self.tu = fixture_taxonomy()

    def test_hand_calculated_metrics_and_cache_safety(self):
        counts = qs.aggregate_lca_mapping(self.mapping)
        original = dict(counts)
        expected = {
            "total_kmers": 17, "unclassified_kmers": 5,
            "classified_kmers": 10, "masked_kmers": 2,
            "total_distance_moved": 20, "average_distance_moved": 2,
            "topology_scale": 2, "total_normalized_distance": 10,
            "average_normalized_distance": 1, "total_upward_loss": 11,
            "average_upward_loss": 1.1, "total_lateral_shift": 9,
            "average_lateral_shift": 0.9, "distance_variance": 3,
        }
        scales = {}
        qs.update_topology_scales(self.tu, pd.Series([13]), scales)
        self.assertEqual(scales, {13: 2})
        for size in (0, 1, 100):
            cache = OrderedDict()
            for _ in range(2):
                self.assertEqual(qs.distance_metrics(self.tu, 13, counts, scales[13], cache, size), expected)
                self.assertEqual(counts, original)

    def test_parsed_ancestor_hits_add_weighted_upward_distance(self):
        # Label 13 -> genus 11 -> family 10: two edges upward per family hit.
        frame = self.frame(["a"], ["10:4 13:3 10:6 0:5 A:2"])
        normalized = qs.normalize_chunk(frame, fill_state(), 2)
        counts = normalized.iloc[0]["lca_counts"]
        self.assertEqual(counts, {10: 10, 13: 3, 0: 5, -1: 2})
        metrics = qs.distance_metrics(self.tu, 13, counts, 2)
        self.assertEqual(metrics["total_upward_loss"], 20)
        self.assertEqual(metrics["total_lateral_shift"], 0)
        self.assertEqual(metrics["total_distance_moved"], 20)
        self.assertAlmostEqual(metrics["average_upward_loss"], 20 / 13)
        self.assertAlmostEqual(metrics["average_distance_moved"], 20 / 13)
        self.assertEqual(metrics["total_normalized_distance"], 10)
        self.assertAlmostEqual(metrics["average_normalized_distance"], 10 / 13)

    def test_each_direction_and_zero_classified(self):
        for taxon, upward, lateral in [(13, 0, 0), (11, 1, 0), (16, 0, 1), (14, 1, 1), (15, 2, 2), (1, 3, 0)]:
            with self.subTest(taxon=taxon):
                result = qs.distance_metrics(self.tu, 13, {taxon: 4}, 2)
                self.assertEqual(result["average_upward_loss"], upward)
                self.assertEqual(result["average_lateral_shift"], lateral)
                self.assertEqual(result["average_distance_moved"], upward + lateral)
                self.assertEqual(result["distance_variance"], 0)
        for counts in ({}, {0: 7, -1: 3}):
            result = qs.distance_metrics(self.tu, 13, counts, 2)
            self.assertEqual(result["classified_kmers"], 0)
            for name, value in result.items():
                if name not in {"total_kmers", "unclassified_kmers", "masked_kmers", "topology_scale"}:
                    self.assertEqual(value, 0)

    def frame(self, accessions, mappings, index=None):
        return pd.DataFrame({
            "accession": accessions, "label_taxon": [13] * len(accessions),
            "a2t_taxon": [13] * len(accessions), "lca_mapping": mappings,
        }, index=index)

    def test_fill_across_chunks_eviction_and_disabled_cache(self):
        for size in (0, 1, 2):
            with self.subTest(size=size):
                state = fill_state()
                first = self.frame(["a", "b"], [self.mapping, "15:2"], [5, 9])
                qs.normalize_chunk(first, state, size)
                second = self.frame(["b", "a", "missing"], [pd.NA] * 3, [12, 15, 20])
                result = qs.normalize_chunk(second, state, size)
                self.assertEqual(result.index.tolist(), [12, 15] if size == 2 else [12])
                self.assertEqual(result.loc[12, "lca_counts"], {15: 2})
                self.assertNotIn("lca_mapping", result)
                self.assertTrue(all(isinstance(v, dict) for v in state["cache"].values()))
                if size == 2:
                    self.assertEqual(result.loc[15, "lca_counts"], qs.aggregate_lca_mapping(self.mapping))
                self.assertIn("lca_mapping", first)  # Caller data is unchanged.

    def test_new_mapping_overrides_old_counts_and_missing_taxa_drop(self):
        frame = self.frame(["a"] * 4, ["13:2", "15:3", pd.NA, "13:1"])
        frame.loc[3, "label_taxon"] = np.nan
        result = qs.normalize_chunk(frame, fill_state(), 2)
        self.assertEqual(result["lca_counts"].tolist(), [{13: 2}, {15: 3}, {15: 3}])
        empty = qs.normalize_chunk(frame.iloc[:0], fill_state(), 2)
        self.assertTrue(empty.empty)
        self.assertNotIn("lca_mapping", empty)

    def test_long_mapping_does_not_infer_fixed_width_unicode(self):
        # The original raw-list assignment would request about 160 GiB here.
        long_mapping = "13:1 " * 171817 + "13:1000"
        self.assertEqual(len(long_mapping), 859092)
        frame = self.frame(["a"] * 50000, [long_mapping] + ["13:1"] * 49999)
        real_asarray = np.asarray

        def guarded_asarray(value, *args, **kwargs):
            if isinstance(value, list) and len(value) == 50000 and all(isinstance(v, str) for v in value):
                dtype = kwargs.get("dtype", args[0] if args else None)
                if dtype is None or np.dtype(dtype).kind == "U":
                    raise AssertionError("Unsafe fixed-width string conversion")
            return real_asarray(value, *args, **kwargs)

        with patch.object(np, "asarray", guarded_asarray):
            result = qs.normalize_chunk(frame, fill_state(), 0)
        self.assertEqual(result.loc[0, "lca_counts"], {13: 172817})
        self.assertEqual(result.loc[49999, "lca_counts"], {13: 1})
        self.assertNotIn("lca_mapping", result)

    def test_cli_five_and_six_columns_preserve_all_scores(self):
        five = f"C\tNC_000001.1\t15\t50\t{self.mapping}\n"
        six = f"C\tNC_000001.1\t13\t50\t{self.mapping}\t13\nC\tNC_000001.1\t15\t50\t\t13\n"
        # Invalid accessions exercise the nonconsecutive index after filtering.
        invalid = "U\tnot_an_accession\t0\t50\t0:16\n"
        with tempfile.TemporaryDirectory() as directory:
            for data, labels in [(invalid + five, [13]), (six, [13, 15])]:
                for chunksize in (1, 10):
                    with self.subTest(labels=labels, chunksize=chunksize):
                        source = Path(directory) / "input.tsv"
                        target = Path(directory) / "scores.csv"
                        source.write_text(data)
                        args = argparse.Namespace(kraken_output=str(source), results=str(target),
                            chunksize=chunksize, verbose=False, lca_cache_size=20000,
                            lca_fill_cache_size=2, movement_cache_size=2)
                        with patch.object(qs, "parse_args", return_value=args), patch.object(qs, "taxutils", return_value=self.tu), patch.object(self.tu, "load_a2t"):
                            qs.main()
                        actual = pd.read_csv(target)
                        self.assertEqual(actual.columns.tolist(), qs.OUTPUT_COLUMNS)
                        self.assertEqual(actual["label_taxon"].tolist(), labels)
                        self.assertEqual(actual["a2t_taxon"].tolist(), [13] * len(labels))
                        for row, label in zip(actual.to_dict("records"), labels):
                            expected = qs.distance_metrics(self.tu, label, qs.aggregate_lca_mapping(self.mapping), 2)
                            for metric, value in expected.items():
                                self.assertAlmostEqual(row[metric], value)

    def test_legacy_cache_option_still_accepted(self):
        with patch("sys.argv", ["quality_score.py", "-i", "in", "-o", "out", "--lca-cache-size", "0"]):
            self.assertEqual(qs.parse_args().lca_cache_size, 0)


if __name__ == "__main__":
    unittest.main()
