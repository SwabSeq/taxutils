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
        (self.root / "targets.json").write_text('{"pathogens":{"test":13}}')
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
        self.assertEqual(first.target_taxa, second.target_taxa)
        self.assertEqual(first.get_branch(13), second.get_branch(13))
        self.assertEqual(first.get_subtree(2), second.get_subtree(2))
        self.assertEqual(first.get_ancestor(13, "D"), second.get_ancestor(13, "D"))

    def test_compressed_object_reuse_skips_all_taxonomy_setup(self):
        raw = {p: (p.read_bytes(), p.stat().st_mtime_ns) for p in (
            self.source, self.root / "names.dmp", self.root / "targets.json",
        )}
        with patch.dict(os.environ, {"TAXUTILS_GLOBALS": str(self.root)}):
            first = constructor.taxutils(low_memory=False)
        self.assertEqual(self.cache.read_bytes()[:2], b"\x1f\x8b")
        with (
            patch.object(constructor, "build_names", side_effect=AssertionError("names rebuilt")),
            patch.object(constructor, "build_nodes", side_effect=AssertionError("nodes rebuilt")),
            patch.object(constructor, "build_parent", side_effect=AssertionError("parents rebuilt")),
            patch.object(constructor, "build_target_taxa", side_effect=AssertionError("targets rebuilt")),
            patch.object(constructor.TaxonomicUtils, "__post_init__", side_effect=AssertionError("setup repeated")),
        ):
            second = self.load()
        self.assert_same_taxonomy(first, second)
        for path, expected in raw.items():
            self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), expected)
        first.names[13] = "Changed by caller"
        first.nodes.loc[first.nodes["taxon"] == 13, "rank_code"] = "F"
        first.target_taxa.clear()
        third = self.load()
        self.assert_same_taxonomy(second, third)

    def test_incompatible_corrupt_and_truncated_caches_rebuild(self):
        expected = self.load()
        installed = self.cache.read_bytes()
        with gzip.open(self.cache, "rb") as cached:
            header = pickle.load(cached)
            tu = pickle.load(cached)
        invalid = [b"broken gzip", installed[:-8]]
        invalid.append(gzip.compress(pickle.dumps((-1, header[1])) + pickle.dumps(tu)))
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

    def test_changed_nodes_names_and_targets_invalidate(self):
        self.load()
        self.source.write_text(self.source.read_text() + "16\t|\t2\t|\tspecies\t|\n")
        self.assertIn(16, self.load().nodes["taxon"].values)
        names_path = self.root / "names.dmp"
        names_path.write_text(names_path.read_text().replace("Alpha", "New alpha"))
        self.assertEqual(self.load().names[13], "New alpha")
        (self.root / "targets.json").write_text('{"pathogens":{"test":15}}')
        self.assertIn(15, self.load().target_taxa)
        stat = self.source.stat()
        os.utime(self.source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000))
        with patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build:
            self.load()
        build.assert_called_once()

    def test_custom_target_path_and_content(self):
        self.load()
        targets = self.root / "custom.json"
        targets.write_text('{"pathogens":{"test":15}}')
        custom = self.load(targets_json=targets)
        self.assertIn(15, custom.target_taxa)
        self.assertNotIn(13, custom.target_taxa)
        with patch.object(constructor, "build_target_taxa", side_effect=AssertionError("rebuilt")):
            self.assert_same_taxonomy(custom, self.load(targets_json=targets))
        targets.write_text('{"pathogens":{"test":13}}')
        self.assertIn(13, self.load(targets_json=targets).target_taxa)
        self.assertIn(13, self.load().target_taxa)

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
            patch.object(constructor, "download_targets"),
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
            patch.object(constructor, "download_targets"),
            patch.object(constructor, "build_nodes", wraps=constructor.build_nodes) as build,
        ):
            self.load(refresh=True)
        download.assert_called_once()
        build.assert_called_once()


if __name__ == "__main__":
    unittest.main()
