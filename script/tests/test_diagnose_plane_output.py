"""Synthetic peaks, time histories, ownership, and malformed input checks."""

import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import diagnose_plane_output as diagnose


class DiagnosePlaneTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.file = self.root / "planedisplacements.0"
        self.out = self.root / "diagnostic"
        self.a = np.zeros((5, 4, 3, 3))
        self.a[:, 0, 1, 0] = np.arange(5) ** 2
        self.a[2, 1, 1, 0] = 0.5
        self.a[1, 3, 2, 1] = -9
        self.a[3, 2, 1, 2] = 0.7
        self.a.tofile(self.file)

    def run_cli(self, *args):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return diagnose.main([str(self.file), "--nx", "4", "--ny", "3", "--dt", "0.1",
                                  "--expected-snapshots", "5", "--top", "1",
                                  "--output-dir", str(self.out), *args])

    def rows(self, name):
        with (self.out / name).open() as stream:
            return list(csv.DictReader(stream))

    def test_actual_peaks_signed_values_and_complete_histories(self):
        self.assertEqual(self.run_cli("--points", "1,0", "--first-time", "2"), 0)
        rows = self.rows("peak_points.csv")
        x, y, z = rows
        self.assertEqual((x["i"], x["j"], x["snapshot"]), ("0", "1", "4"))
        self.assertAlmostEqual(float(x["time_s"]), 2.4)
        self.assertEqual(float(y["value_m"]), -9)
        self.assertEqual(float(z["value_m"]), 0.7)
        self.assertEqual(x["medium"], "unknown")
        self.assertEqual(x["cpml"], "-1")
        histories = self.rows("histories.csv")
        for row in histories:
            step, i, j = (int(row[k]) for k in ("snapshot", "i", "j"))
            np.testing.assert_array_equal([float(row[k]) for k in ("ux_m", "uy_m", "uz_m")],
                                          self.a[step, i, j])
        self.assertEqual(len(histories), 5 * 4)

    def test_boundary_band_is_separate_from_exact_cpml_membership(self):
        meta = self.root / "plane_point_metadata.0.csv"
        meta.write_text("plane,i,j,medium,cpml,x_m,y_m,z_m,valid,acoustic_reconstruction,pressure_defined\n"
                        "0,0,1,acoustic,1,0,10,30,1,1,1\n"
                        "0,1,1,solid,1,10,10,20,1,0,0\n"
                        "0,2,1,solid,0,20,10,20,1,0,0\n")
        self.assertEqual(self.run_cli("--spacing", "10", "10", "--domain", "30", "20",
                                      "--cpml-thickness", "5", "--metadata", str(meta)), 0)
        summary = json.loads((self.out / "summary.json").read_text())
        groups = summary["groups"]
        self.assertEqual(groups["lateral_boundary_band"]["points"], 10)
        self.assertEqual(groups["geometric_interior"]["peak_abs_m"], [0.5, 0, 0.7])
        self.assertEqual(groups["cpml"]["points"], 2)
        self.assertEqual(groups["unknown"]["points"], 9)
        x = self.rows("peak_points.csv")[0]
        self.assertEqual((x["medium"], x["cpml"]), ("acoustic", "1"))

    def test_equal_peaks_use_first_time_and_lowest_point_index(self):
        self.a.fill(1)
        self.a.tofile(self.file)
        self.assertEqual(self.run_cli(), 0)
        for row in self.rows("peak_points.csv"):
            self.assertEqual((row["i"], row["j"], row["snapshot"]), ("0", "0", "0"))

    def test_rejects_partial_truncated_and_nonfinite_inputs(self):
        self.file.write_bytes(b"partial")
        self.assertEqual(self.run_cli(), 2)
        self.a[:-1].tofile(self.file)
        self.assertEqual(self.run_cli(), 2)
        for bad in (np.nan, np.inf):
            self.a[-1, 1, 2, 2] = bad
            self.a.tofile(self.file)
            self.assertEqual(self.run_cli(), 2)
        self.assertFalse(self.out.exists())

    def test_big_endian_and_existing_output_directory(self):
        self.a.astype(">f8").tofile(self.file)
        self.assertEqual(self.run_cli("--endian", "big"), 0)
        self.assertEqual(float(self.rows("peak_points.csv")[0]["value_m"]), 16)
        self.assertEqual(self.run_cli("--endian", "big"), 2)

    def test_duplicate_metadata_and_geometry_mismatch(self):
        meta = self.root / "metadata.csv"
        header = "plane,i,j,medium,cpml,x_m,y_m,z_m,valid,acoustic_reconstruction,pressure_defined\n"
        row = "0,0,1,acoustic,1,1,10,30,1,1,1\n"
        meta.write_text(header + row + row)
        self.assertEqual(self.run_cli("--metadata", str(meta)), 2)
        meta.write_text(header + row)
        self.assertEqual(self.run_cli("--metadata", str(meta), "--spacing", "10", "10",
                                      "--domain", "30", "20", "--cpml-thickness", "5"), 2)


if __name__ == "__main__":
    unittest.main()
