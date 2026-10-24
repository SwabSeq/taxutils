import gzip
import importlib
import os
import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from taxutils import cache

constructor = importlib.import_module("taxutils.taxutils")


class TaxonomyCacheTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "nodes.dmp"
        self.cache = self.root / "taxutils.pkl.gz"
        self.source.write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "13\t|\t2\t|\tspecies\t|\n"
            "15\t|\t2\t|\tspecies\t|\n"
        )
        (self.root / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
            "13\t|\tAlpha\t|\t\t|\tscientific name\t|\n"
            "15\t|\tBeta\t|\t\t|\tscientific name\t|\n"
        )
        for name in ("ensure_a2t_db", "prepare_alternative_mappings"):
            mock = patch.object(constructor, name)
            mock.start()
            self.addCleanup(mock.stop)

    def load(self, **kwargs):
        return constructor.taxutils(save_folder=self.root, low_memory=False, **kwargs)

    def assert_same_taxonomy(self, first, second):
        pd.testing.assert_frame_equal(first.nodes, second.nodes)
        self.assertEqual(first.names, second.names)
        self.assertEqual(first.parent, second.parent)
        self.assertEqual(first.get_branch(13), second.get_branch(13))
        self.assertEqual(first.get_subtree(2), second.get_subtree(2))
        self.assertEqual(first.get_ancestor(13, "D"), second.get_ancestor(13, "D"))

    def test_compressed_object_reuse_skips_all_taxonomy_setup(self):
        raw = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (
            self.source, self.root / "names.dmp",
        )}
        with patch.dict(os.environ, {"TAXUTILS_GLOBALS": str(self.root)}):
            first = constructor.taxutils(low_memory=False)
        self.assertEqual(self.cache.read_bytes()[:2], b"\x1f\x8b")
        with (
            patch.object(constructor, "build_names", side_effect=AssertionError("names rebuilt")),
            patch.object(constructor, "build_nodes", side_effect=AssertionError("nodes rebuilt")),
            patch.object(constructor, "build_parent", side_effect=AssertionError("parents rebuilt")),
            patch.object(constructor.TaxonomicUtils, "__post_init__", side_effect=AssertionError("setup repeated")),
        ):
            second = self.load()
        self.assert_same_taxonomy(first, second)
        for path, expected in raw.items():
            self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), expected)
        first.names[13] = "Changed by caller"
        first.nodes.loc[first.nodes["taxon"] == 13, "rank_code"] = "F"
        third = self.load()
        self.assert_same_taxonomy(second, third)

    def test_categorical_ranks_preserve_values_and_comparisons(self):
        for low_memory in (True, False):
            with self.subTest(low_memory=low_memory):
                tu = constructor.taxutils(save_folder=self.root, low_memory=low_memory)
                nodes = tu.nodes
                for column in ("rank", "rank_code", "rank_base", "new_rank"):
                    self.assertIsInstance(nodes[column].dtype, pd.CategoricalDtype)
                    self.assertFalse(nodes[column].cat.ordered)
                    self.assertFalse(nodes[column].isna().any())
                for column in ("taxon", "parent", "rank_idx"):
                    self.assertTrue(pd.api.types.is_integer_dtype(nodes[column]))
                self.assertEqual(nodes["rank"].tolist(), ["no rank", "superkingdom", "species", "species"])
                self.assertEqual(nodes["rank_code"].tolist(), ["R", "D", "S", "S"])
                self.assertEqual(nodes["new_rank"].tolist(), ["root", "domain", "species", "species"])
                self.assertEqual((nodes["rank"] != nodes["new_rank"]).tolist(), [True, True, False, False])
                self.assertTrue((nodes["rank_code"] == nodes["rank_base"]).all())
                self.assertEqual(nodes["rank"].value_counts()["species"], 2)

    def test_object_rank_cache_rebuilds_with_categoricals(self):
        self.load()
        with gzip.open(self.cache, "rb") as cached:
            header = pickle.load(cached)
            tu = pickle.load(cached)
        tu.nodes = tu.nodes.astype({column: object for column in (
            "rank", "rank_code", "rank_base", "new_rank",
        )})
        self.cache.write_bytes(gzip.compress(
            pickle.dumps((1, header)) + pickle.dumps(tu)
        ))
        with patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build:
            rebuilt = self.load()
        build.assert_called_once()
        self.assertIsInstance(rebuilt.nodes["rank"].dtype, pd.CategoricalDtype)
        self.assert_same_taxonomy(rebuilt, self.load())

    def test_legacy_header_reuses_compatible_taxonomy(self):
        expected = self.load()
        with gzip.open(self.cache, "rb") as cached:
            sources = pickle.load(cached)
            tu = pickle.load(cached)
        self.cache.write_bytes(gzip.compress(
            pickle.dumps((2, sources)) + pickle.dumps(tu)
        ))
        with patch.object(constructor, "build_nodes", side_effect=AssertionError("nodes rebuilt")):
            self.assert_same_taxonomy(expected, self.load())

    def test_incompatible_corrupt_and_truncated_caches_rebuild(self):
        expected = self.load()
        installed = self.cache.read_bytes()
        with gzip.open(self.cache, "rb") as cached:
            header = pickle.load(cached)
            tu = pickle.load(cached)
        invalid = [b"broken gzip", installed[:-8]]
        invalid.append(gzip.compress(pickle.dumps(header[:-1]) + pickle.dumps(tu)))
        invalid.append(gzip.compress(pickle.dumps(header) + pickle.dumps(None)))
        tu.nodes = tu.nodes.drop(columns="rank_code")
        invalid.append(gzip.compress(pickle.dumps(header) + pickle.dumps(tu)))
        for content in invalid:
            with self.subTest(content=content[:20]):
                self.cache.write_bytes(content)
                with patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build:
                    actual = self.load()
                build.assert_called_once()
                self.assert_same_taxonomy(expected, actual)

    def test_changed_nodes_and_names_invalidate(self):
        self.load()
        self.source.write_text(self.source.read_text() + "16\t|\t2\t|\tspecies\t|\n")
        self.assertIn(16, self.load().nodes["taxon"].values)
        names_path = self.root / "names.dmp"
        names_path.write_text(names_path.read_text().replace("Alpha", "New alpha"))
        self.assertEqual(self.load().names[13], "New alpha")
        stat = self.source.stat()
        os.utime(self.source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
        with patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build:
            self.load()
        build.assert_called_once()

    def test_target_free_object_and_legacy_cache(self):
        tu = self.load()
        self.assertFalse(hasattr(tu, "target_taxa"))
        self.assertFalse((self.root / "targets.json").exists())
        tu.target_taxa = [13]
        sources = cache.taxonomy_sources(self.root / "names.dmp", self.source)
        cache.save_taxonomy_cache(self.cache, sources, tu)
        with patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build:
            restored = self.load()
        build.assert_called_once()
        self.assertFalse(hasattr(restored, "target_taxa"))

    def test_current_options_accessions_and_backend_preparation(self):
        with patch.object(constructor, "build_a2t", return_value={"NC_000001.1": 13}):
            first = self.load(accessions=["NC_000001.1"], threads=1)
        self.assertEqual(first.a2t["NC_000001.1"], 13)
        with (
            patch.object(constructor, "ensure_a2t_db") as ensure,
            patch.object(constructor, "prepare_alternative_mappings") as alternatives,
            patch.object(constructor, "build_a2t", return_value={"NC_000002.1": 15}) as lookup,
        ):
            second = self.load(accessions=["NC_000002.1"], canonical=False, wgs=True, threads=2)
        self.assertEqual(second.a2t, {"NC_000002.1": 15, 0: "unclassified"})
        self.assertFalse(second._canonical)
        self.assertTrue(second._wgs)
        self.assertEqual(second._threads, 2)
        ensure.assert_called_once_with(save_folder=str(self.root), wgs=True, refresh=False, threads=2)
        alternatives.assert_called_once_with(
            save_folder=str(self.root), canonical=False, low_memory=False,
            wgs=True, refresh=False, threads=2,
        )
        lookup.assert_called_once_with(
            ["NC_000002.1"], save_folder=str(self.root), low_memory=False,
            wgs=True, threads=2, canonical=False,
        )
        third = self.load()
        self.assertIsNone(third.a2t)
        self.assertTrue(third._canonical)
        self.assertFalse(third._wgs)
        self.assertIsNone(third._threads)

    def test_failed_cache_install_preserves_installed_object(self):
        first = self.load()
        installed = self.cache.read_bytes()
        with (
            patch.object(constructor, "download_taxdump"),
            patch.object(cache.os, "replace", side_effect=OSError("write failed")),
        ):
            second = self.load(refresh=True)
        self.assert_same_taxonomy(first, second)
        self.assertEqual(self.cache.read_bytes(), installed)
        self.assertEqual(list(self.root.glob(".taxutils-*.tmp")), [])

    def test_low_memory_ignores_cache_and_legacy_nodes_cache(self):
        self.load()
        installed = self.cache.read_bytes()
        (self.root / "cnodes.dmp").write_bytes(b"obsolete")
        with patch.object(constructor, "load_taxonomy_cache", side_effect=AssertionError("cache read")):
            constructor.taxutils(save_folder=self.root, low_memory=True)
            self.assertEqual(self.cache.read_bytes(), installed)
            self.cache.unlink()
            constructor.taxutils(save_folder=self.root, low_memory=True)
            self.assertFalse(self.cache.exists())
        self.load()
        self.assertTrue(self.cache.exists())
        self.assertEqual((self.root / "cnodes.dmp").read_bytes(), b"obsolete")

    def test_refresh_forces_rebuild(self):
        self.load()
        with (
            patch.object(constructor, "download_taxdump") as download,
            patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build,
        ):
            self.load(refresh=True)
        download.assert_called_once()
        build.assert_called_once()


if __name__ == "__main__":
    unittest.main()
