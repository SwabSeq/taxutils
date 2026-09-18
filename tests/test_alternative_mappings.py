"""Public Python integration tests for native viral overrides."""
import gzip
from pathlib import Path
import sqlite3
import tempfile
import unittest

from taxutils import taxutils
from taxutils.backend import TaxutilsBackendError


class AlternativeMappingTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        (self.root / "nodes.dmp").write_text(
            "1 | 1 | no rank |\n10 | 1 | species |\n"
            "11 | 10 | no rank |\n12 | 11 | no rank |\n13 | 12 | no rank |\n"
        )
        (self.root / "names.dmp").write_text(
            "1 | root | | scientific name |\n10 | Influenza A virus | | scientific name |\n"
            "11 | intermediate | | scientific name |\n12 | subtype | | scientific name |\n"
            "13 | Influenza A virus (A/Test/1/2020(H1N1)) | | scientific name |\n"
        )
        (self.root / "targets.json").write_text('{"pathogens": {}}')
        self.write_gzip("nucl_gb.accession2taxid.gz",
                        "accession\taccession.version\ttaxid\tgi\n"
                        "NC_000001\tNC_000001.1\t10\t0\n"
                        "NC_000002\tNC_000002.1\t10\t0\n"
                        "NC_000003\tNC_000003.1\t10\t0\n")
        self.write_gzip("viral.metadata.csv.gz",
                        "#Accession,Genotype,Strain\n"
                        "NC_000001.1,,A/Test/1/2020\n"
                        "NC_000002.1,H1N1,unknown\n")

    def write_gzip(self, name, content):
        with gzip.open(self.root / name, "wt") as out:
            out.write(content)

    def test_constructor_extend_modes_and_reverse(self):
        for low_memory in (True, False):
            with self.subTest(low_memory=low_memory):
                tu = taxutils(save_folder=self.root, low_memory=low_memory,
                              canonical=False, threads=2, accessions=["NC_000001.1"])
                self.assertEqual(tu.a2t["NC_000001.1"], 13)
                self.assertNotIn("strain", tu.nodes.columns)
                self.assertNotIn("subtype", tu.nodes.columns)
                tu.load_a2t([">NC_000002.1 description", "NC_000003.1"], extend=True)
                self.assertEqual({a: tu.a2t[a] for a in ("NC_000001.1", "NC_000002.1", "NC_000003.1")},
                                 {"NC_000001.1": 13, "NC_000002.1": 12, "NC_000003.1": 10})
                self.assertEqual(tu.get_t2a([10]), {"NC_000003.1"})
                self.assertEqual(tu.get_t2a([12, 13]), {"NC_000001.1", "NC_000002.1"})
                tu.load_a2t(["NC_000001.2", "NC_000001.1"], low_memory=not low_memory)
                self.assertEqual(tu.a2t, {"NC_000001.1": 13})
                tu.load_a2t([])
                self.assertEqual(tu.a2t, {})

    def test_canonical_coexists_and_needs_no_metadata(self):
        alt = taxutils(save_folder=self.root, canonical=False, low_memory=False, threads=1)
        alt.load_a2t(["NC_000001.1"])
        (self.root / "viral.metadata.csv.gz").unlink()
        for low_memory in (True, False):
            default = taxutils(save_folder=self.root, low_memory=low_memory, threads=1)
            default.load_a2t(["NC_000001.1"])
            self.assertEqual(default.a2t, {"NC_000001.1": 10})
            self.assertEqual(default.get_t2a([10]), {"NC_000001.1", "NC_000002.1", "NC_000003.1"})
        self.assertFalse((self.root / "viral.metadata.csv.gz").exists())

    def test_metadata_replacement_and_corruption(self):
        tu = taxutils(save_folder=self.root, canonical=False, low_memory=False, threads=1)
        self.write_gzip("viral.metadata.csv.gz", "#Accession,Genotype,Strain\nNC_000001.1,H1N1,unknown\n")
        tu.load_a2t(["NC_000001.1", "NC_000002.1"])
        self.assertEqual(tu.a2t, {"NC_000001.1": 12, "NC_000002.1": 10})
        with sqlite3.connect(self.root / "nucl.accession2taxid.db") as db:
            self.assertEqual(db.execute("SELECT count(*) FROM a2t_overrides").fetchone()[0], 1)
        self.write_gzip("viral.metadata.csv.gz", "broken,header\na,b\n")
        with self.assertRaisesRegex(TaxutilsBackendError, "metadata missing"):
            tu.load_a2t(["NC_000001.1"])


if __name__ == "__main__":
    unittest.main()
