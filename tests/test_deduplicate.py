from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from taxutils.deduplicate_fasta import deduplicate_fasta


class DeduplicateTests(unittest.TestCase):
    def test_cli_output_and_in_place(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "input.fa"
            original = b">NC_000001.1 first\r\nAC\r\n>kraken:taxid|13|NC_000001.1 copy\nTT\n>NC_000001.2\nAC"
            expected = b">NC_000001.1 first\r\nAC\r\n>NC_000001.2\nAC"
            source.write_bytes(original)
            for extra in [["-o", "out.fa"], []]:
                result = subprocess.run(
                    [sys.executable, "-m", "taxutils.cli", "deduplicate",
                     "-i", "input.fa", "--threads", "2", *extra],
                    cwd=root, capture_output=True, text=True,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("kept=2 removed=1", result.stdout)
                self.assertEqual((root / "out.fa").read_bytes(), expected)
                self.assertEqual(source.read_bytes(), original if extra else expected)

    def test_bad_header_leaves_input_and_output_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "input.fa"
            output = Path(directory) / "output.fa"
            original = b">NC_000001.1\nAA\n>unknown\nTT\n"
            source.write_bytes(original)
            output.write_bytes(b"existing")
            for destination in [None, output]:
                with self.assertRaisesRegex(ValueError, "No accession found"):
                    deduplicate_fasta(source, destination)
                self.assertEqual(source.read_bytes(), original)
                self.assertEqual(output.read_bytes(), b"existing")

    def test_requires_input(self):
        result = subprocess.run(
            [sys.executable, "-m", "taxutils.cli", "deduplicate", "-o", "unused.fa"],
            capture_output=True, text=True,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--input", result.stderr)


if __name__ == "__main__":
    unittest.main()
