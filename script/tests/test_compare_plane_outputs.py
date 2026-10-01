"""Synthetic binary fixtures for the plane-output equivalence checker."""

import contextlib
import io
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import compare_plane_outputs as compare


class ComparePlanesTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.ref = Path(self.temp.name) / "reference"
        self.test = Path(self.temp.name) / "candidate"
        self.ref.mkdir()
        self.test.mkdir()
        self.a = np.arange(4 * 2 * 3 * 3, dtype=np.float64).reshape(4, 2, 3, 3) / 100

    def write_pair(self, a=None, b=None, index=0):
        a = self.a if a is None else a
        b = a if b is None else b
        ref = self.ref / f"planedisplacements.{index}"
        test = self.test / f"planedisplacements.{index}"
        a.tofile(ref)
        b.tofile(test)
        return ref, test

    def cli(self, *options):
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            result = compare.main([str(self.ref), str(self.test),
                                   "--shapes", "2x3", *options])
        return result, output.getvalue(), errors.getvalue()

    def test_identical_files_and_multiple_chunks(self):
        self.write_pair()
        code, out, _ = self.cli("--expected-snapshots", "4", "--chunk-points", "5")
        self.assertEqual(code, 0)
        self.assertIn("byte-identical=yes", out)
        self.assertIn("exact finite values=yes", out)
        self.assertIn("OVERALL: PASS", out)

    def test_difference_location_across_frame_and_chunk_boundaries(self):
        b = self.a.copy()
        b[2, 1, 2, 1] += 0.25
        ref, test = self.write_pair(b=b)
        identical, stats = compare.compare_plane(ref, test, 2, 3, 4, chunk_points=5)
        self.assertFalse(identical)
        self.assertEqual(stats[1].max_location[:3], (2, 1, 2))
        self.assertEqual(stats[1].outside_tolerance, 1)
        self.assertAlmostEqual(stats[1].max_abs_error, 0.25)
        self.assertTrue(stats[0].passed and stats[2].passed)
        self.assertEqual(self.cli()[0], 1)

    def test_tolerance_uses_each_reference_value_including_zero(self):
        a = np.zeros_like(self.a)
        a[0, 0, 0, 0] = 1000
        b = a.copy()
        b[0, 0, 0, 0] += 0.0009  # accepted by rtol
        b[3, 1, 2, 0] = 2e-9   # zero reference: rejected by atol
        ref, test = self.write_pair(a, b)
        _, stats = compare.compare_plane(ref, test, 2, 3, 4,
                                         rtol=1e-6, atol=1e-10, chunk_points=5)
        self.assertEqual(stats[0].finite_differences, 2)
        self.assertEqual(stats[0].outside_tolerance, 1)
        self.assertEqual(stats[0].first_failure[:3], (3, 1, 2))
        code, out, _ = self.cli("--rtol", "1e-6", "--atol", "3e-9")
        self.assertEqual(code, 0)
        self.assertIn("exact finite values=no", out)

    def test_matching_nan_or_infinity_never_pass(self):
        for value in (np.nan, np.inf, -np.inf):
            with self.subTest(value=value):
                a = self.a.copy()
                a[-1, -1, -1, 2] = value
                ref, test = self.write_pair(a)
                identical, stats = compare.compare_plane(ref, test, 2, 3, 4)
                self.assertTrue(identical)
                self.assertFalse(stats[2].passed)
                self.assertEqual(stats[2].nonfinite_reference, 1)
                self.assertEqual(stats[2].nonfinite_candidate, 1)
                self.assertEqual(self.cli()[0], 1)

    def test_signed_zero_is_numerically_equal_but_not_byte_identical(self):
        a = np.zeros_like(self.a)
        b = a.copy()
        b[0, 0, 0, 0] = -0.0
        self.write_pair(a, b)
        code, out, _ = self.cli()
        self.assertEqual(code, 0)
        self.assertIn("byte-identical=no", out)
        self.assertIn("exact finite values=yes", out)

    def test_big_endian(self):
        b = self.a.copy()
        b[1, 0, 1, 0] += 0.125
        ref, test = self.write_pair(self.a.astype(">f8"), b.astype(">f8"))
        _, stats = compare.compare_plane(ref, test, 2, 3, 4, endian="big")
        self.assertAlmostEqual(stats[0].max_abs_error, 0.125)
        self.assertEqual(stats[0].max_location[:3], (1, 0, 1))

    def test_truncated_snapshot_and_mismatched_counts(self):
        for candidate in (self.a.ravel()[:-1], self.a[:-1]):
            with self.subTest(size=candidate.size):
                self.write_pair(b=candidate)
                self.assertEqual(self.cli()[0], 2)

    def test_expected_count_rejects_equally_truncated_files(self):
        self.write_pair(self.a[:-1])
        code, _, error = self.cli("--expected-snapshots", "4")
        self.assertEqual(code, 2)
        self.assertIn("3 snapshots; expected 4", error)

    def test_missing_and_extra_planes(self):
        ref, test = self.write_pair()
        test.unlink()
        self.assertEqual(self.cli()[0], 2)
        self.write_pair()
        self.write_pair(index=1)
        self.assertEqual(self.cli()[0], 2)

    def test_zero_length_file(self):
        self.write_pair(np.array([], dtype=np.float64))
        self.assertEqual(self.cli()[0], 2)

    def test_all_planes_are_checked_after_numerical_failure(self):
        b = self.a.copy()
        b[0, 0, 0, 0] += 1
        self.write_pair(b=b)
        self.write_pair(index=1)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            code = compare.main([str(self.ref), str(self.test),
                                 "--shapes", "2x3", "2x3"])
        self.assertEqual(code, 1)
        self.assertIn("planedisplacements.1: PASS", out.getvalue())

    def test_difference_overflow_does_not_pass(self):
        a = np.full_like(self.a, 1e308)
        ref, test = self.write_pair(a, -a)
        with np.errstate(all="raise"):
            _, stats = compare.compare_plane(ref, test, 2, 3, 4)
        self.assertTrue(all(not stat.passed for stat in stats))

    def test_file_changed_during_read_is_rejected(self):
        ref, test = self.write_pair()
        original_read = np.fromfile
        changed = False

        def read_and_touch(*args, **kwargs):
            nonlocal changed
            result = original_read(*args, **kwargs)
            if not changed:
                before = ref.stat()
                os.utime(ref, ns=(before.st_atime_ns, before.st_mtime_ns + 1000000000))
                changed = True
            return result

        with mock.patch.object(compare.np, "fromfile", side_effect=read_and_touch):
            with self.assertRaisesRegex(ValueError, "File changed during comparison"):
                compare.compare_plane(ref, test, 2, 3, 4, chunk_points=5)


if __name__ == "__main__":
    unittest.main()
