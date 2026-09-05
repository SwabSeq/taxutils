import gzip
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest


class NativeBackendTests(unittest.TestCase):
    @classmethod
    def write_accession_fixture(cls):
        with gzip.open(cls.root / "nucl_gb.accession2taxid.gz", "wt") as output:
            output.write("accession\taccession.version\ttaxid\tgi\n")
            output.write("NC_000001\tNC_000001.1\t13\t0\n")
            output.write("NC_000002\tNC_000002.1\t15\t0\n")

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        os.environ["TAXUTILS_GLOBALS"] = str(cls.root)

        (cls.root / "names.dmp").write_text(
            "1\t|\troot\t|\t\t|\tscientific name\t|\n"
            "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
            "13\t|\tTestus alpha\t|\t\t|\tscientific name\t|\n"
            "15\t|\tTestus beta\t|\t\t|\tscientific name\t|\n"
        )
        (cls.root / "nodes.dmp").write_text(
            "1\t|\t1\t|\tno rank\t|\n"
            "2\t|\t1\t|\tsuperkingdom\t|\n"
            "13\t|\t2\t|\tspecies\t|\n"
            "15\t|\t2\t|\tspecies\t|\n"
        )
        (cls.root / "targets.json").write_text('{"pathogens":{"example":13}}')
        cls.write_accession_fixture()

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_backend_info_has_stable_shape(self):
        from taxutils import backend_info

        info = backend_info()
        self.assertEqual(info["selected"], "rust")
        self.assertEqual(info["api_version"], 2)
        self.assertEqual(
            set(info),
            {"selected", "rust_version", "api_version"},
        )

    def test_fasta_commands_preserve_outputs_and_stats(self):
        from taxutils.clean_headers import clean_fasta_headers
        from taxutils.extract_headers import extract_accessions
        from taxutils.grep_fasta import grep_fasta

        source = self.root / "input.fa"
        source.write_text(
            ">NC_000001.1 alpha\nAA\n"
            ">NC_000002.1 beta\nCC\n"
        )
        accessions = self.root / "accessions.txt"
        self.assertEqual(extract_accessions(source, accessions, batch_size=1), 2)
        self.assertEqual(accessions.read_text(), "NC_000001.1\nNC_000002.1\n")

        hits = self.root / "hits.fa"
        self.assertEqual(
            grep_fasta(source, "NC_000002.1", hits, batch_size=64),
            {"requested": 1, "scanned": 2, "matched": 1, "missing_accession": 0},
        )
        self.assertEqual(hits.read_text(), ">NC_000002.1 beta\nCC\n")

        cleaned = self.root / "clean.fa"
        self.assertEqual(clean_fasta_headers(source, cleaned), cleaned)
        self.assertEqual(
            cleaned.read_text(),
            ">NC_000001.1\nAA\n>NC_000002.1\nCC\n",
        )

    def test_database_build_can_discard_compressed_source(self):
        from taxutils import taxutils

        database = self.root / "nucl.accession2taxid.db"
        source = self.root / "nucl_gb.accession2taxid.gz"
        database.unlink(missing_ok=True)
        self.write_accession_fixture()
        try:
            tu = taxutils(low_memory=False, keep_accession_downloads=False)
            self.assertTrue(database.exists())
            self.assertFalse(source.exists())

            # Reusing the complete database must not redownload the discarded source.
            tu.load_a2t(["NC_000001.1"])
            self.assertEqual(tu.a2t, {"NC_000001.1": 13})
            self.assertEqual(tu.get_t2a([13]), {"NC_000001.1"})
            self.assertEqual(tu.get_t2a([13, 13, 99999]), {"NC_000001.1"})
            self.assertEqual(tu.get_t2a([13, 15]), {"NC_000001.1", "NC_000002.1"})
            self.assertEqual(tu.get_t2a([99999]), set())
            self.assertFalse(source.exists())

            with sqlite3.connect(database) as connection:
                indexes = {
                    row[0]
                    for row in connection.execute(
                        "SELECT name FROM sqlite_master "
                        "WHERE type = 'index' AND name LIKE 'idx_%'"
                    )
                }
            self.assertEqual(indexes, {"idx_accession", "idx_taxid"})
        finally:
            database.unlink(missing_ok=True)
            self.write_accession_fixture()

    def test_bulk_lookups_and_filter(self):
        from taxutils.filter_fasta import filter_fasta
        from taxutils.taxutils import build_a2t, get_t2a

        self.assertEqual(
            build_a2t(["NC_000001.1", "missing"]),
            {"NC_000001.1": 13},
        )
        self.assertEqual(get_t2a([13]), {"NC_000001.1"})

        source = self.root / "filter.fa"
        output = self.root / "filtered.fa"
        source.write_text(
            ">NC_000001.1 alpha\nAA\n"
            ">NC_000002.1 beta\nCC\n"
            ">unparseable header\nGG\n"
        )
        self.assertEqual(
            filter_fasta(source, output, {13}, filter_mode="keep", batch_size=1),
            {"kept": 1, "removed": 2, "missing_accession": 1, "missing_taxid": 0},
        )
        self.assertEqual(output.read_text(), ">NC_000001.1 alpha\nAA\n")

    def test_missing_extension_is_reported(self):
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import sys; sys.modules['taxutils._rust'] = None; import taxutils",
            ],
            env=os.environ.copy(),
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires its native Rust extension", result.stderr)


if __name__ == "__main__":
    unittest.main()
