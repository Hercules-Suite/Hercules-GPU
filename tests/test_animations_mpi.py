"""Analytical checks for plane layout, velocity, and MPI frame ownership."""

import contextlib
import io
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from script.animations_mpi import (detect_rank_info, load_plane, main, partition,
                                   vertical_velocity)


class PlaneAnimationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.previous_cwd = Path.cwd()
        os.chdir(self.temporary.name)
        Path("planes").mkdir()
        self.plane_file = Path("planes/planedisplacements.0")
        self.options = ["--plane-file", str(self.plane_file), "--nx", "3", "--ny", "2",
                        "--dx", "500", "--dy", "325", "--dt", "0.02"]
        # Write in the C writer's nested order, independently of NumPy reshape.
        # The asymmetric grid and distinct XYZ values detect transposes,
        # mirrors, component-major reads, and accidental use of another component.
        with Path("planes/planedisplacements.0").open("wb") as output:
            for t in (0.0, 0.02, 0.04, 0.06):
                for i in range(3):
                    for j in range(2):
                        uz = 1000 * i + 100 * j + (1 + 10 * i + j) * t
                        output.write(struct.pack("=ddd", -9000 - t, 8000 + t, uz))

    def tearDown(self):
        os.chdir(self.previous_cwd)
        self.temporary.cleanup()

    def load(self, **overrides):
        options = dict(plane_file=self.plane_file, nx=3, ny=2, dt=0.02, dx=500, dy=325)
        options.update(overrides)
        return load_plane(**options)

    def test_binary_point_order_and_plot_coordinates(self):
        plane = self.load()
        data = plane.open()
        self.assertEqual(data.shape, (4, 3, 2, 3))
        np.testing.assert_array_equal(data[0, :, :, 2], [[0, 100], [1000, 1100], [2000, 2100]])
        self.assertEqual(data[0, 2, 1, 0], -9000)
        self.assertEqual(data[0, 2, 1, 1], 8000)
        horizontal, vertical = plane.coordinates()
        np.testing.assert_allclose(horizontal, [0, 0.325])
        np.testing.assert_allclose(vertical, [0, 0.5, 1.0])
        self.assertEqual(plane.dt, 0.02)

    def test_linear_velocity_all_frames_and_no_amplitude_clipping(self):
        plane = self.load()
        expected = [[100, 200], [1100, 1200], [2100, 2200]]
        for frame in range(plane.nt):
            np.testing.assert_allclose(vertical_velocity(plane.open(), frame, plane.dt),
                                       expected, rtol=1e-10)

    def test_central_difference_on_quadratic_displacement(self):
        times = np.arange(5) * 0.02
        data = np.zeros((5, 3, 2, 3))
        data[..., 2] = times[:, None, None] ** 2
        for frame in (1, 2, 3):
            np.testing.assert_allclose(vertical_velocity(data, frame, 0.02),
                                       100 * 2 * times[frame])
        # With only two snapshots, both endpoints use the same secant slope.
        for frame in (0, 1):
            np.testing.assert_allclose(vertical_velocity(data[:2], frame, 0.02), 2.0)

    def test_explicit_file_selects_plane_independently_of_filename(self):
        selected = Path("planes/planedisplacements.7")
        self.plane_file.rename(selected)
        plane = self.load(plane_file=selected, dx=750, dy=150, dt=0.04)
        self.assertEqual((plane.nx, plane.ny, plane.dt), (3, 2, 0.04))
        self.assertEqual(plane.path, selected.resolve())
        horizontal, vertical = plane.coordinates()
        np.testing.assert_allclose(horizontal, [0, 0.15])
        np.testing.assert_allclose(vertical, [0, 0.75, 1.5])
        # --dt is already the saved-snapshot interval: do not multiply by a rate.
        np.testing.assert_allclose(vertical_velocity(plane.open(), 1, plane.dt),
                                   [[50, 100], [550, 600], [1050, 1100]], rtol=1e-10)

    def test_incomplete_bytes_and_insufficient_snapshots(self):
        path = Path("planes/planedisplacements.0")
        original = path.read_bytes()
        for extra in (b"x", b"12345678"):
            path.write_bytes(original + extra)
            with self.assertRaisesRegex(ValueError, "not a multiple"):
                self.load()
        for size in (0, 3 * 2 * 3 * 8):
            path.write_bytes(original[:size])
            with self.assertRaisesRegex(ValueError, "At least two"):
                self.load()

    def test_invalid_sampling_and_nonfinite_velocity(self):
        for name in ("dt", "dx", "dy"):
            for value in (0, -0.001, float("nan"), float("inf")):
                with self.assertRaisesRegex(ValueError, "finite and positive"):
                    self.load(**{name: value})
        for name in ("nx", "ny"):
            for value in (0, 1, -2):
                with self.assertRaisesRegex(ValueError, "at least two"):
                    self.load(**{name: value})
        data = np.zeros((2, 3, 2, 3))
        data[1, 0, 1, 2] = np.nan
        with self.assertRaisesRegex(ValueError, "Non-finite"):
            vertical_velocity(data, 0, 0.02)

    def test_mpi_and_worker_partitions_cover_frames_once(self):
        for frames in (2, 7, 25):
            for ranks in (1, 3, 10):
                owned = []
                for rank in range(ranks):
                    start, stop = partition(frames, rank, ranks)
                    jobs = min(4, stop - start)
                    for job in range(jobs):
                        left, right = partition(stop - start, job, jobs)
                        owned.extend(range(start + left, start + right))
                self.assertEqual(owned, list(range(frames)))
        self.assertEqual(detect_rank_info({"OMPI_COMM_WORLD_RANK": "1",
                                          "OMPI_COMM_WORLD_SIZE": "4",
                                          "SLURM_PROCID": "0", "SLURM_NTASKS": "2"}), (1, 4))
        self.assertEqual(detect_rank_info({"SLURM_PROCID": "2", "SLURM_NTASKS": "3"}), (2, 3))
        self.assertEqual(detect_rank_info({}), (0, 1))
        with self.assertRaisesRegex(ValueError, "world size"):
            detect_rank_info({"PMIX_RANK": "2"})

    def test_check_only_without_parameters_file_and_manual_time_origin(self):
        output = io.StringIO()
        with patch.dict(os.environ, {}, clear=True), contextlib.redirect_stdout(output):
            self.assertFalse(Path("parameters.in").exists())
            main(self.options + ["--check-only"])
            # An unrelated, malformed configuration must never be opened.
            Path("parameters.in").write_text("multiple unrelated planes; do not parse")
            with patch.object(Path, "read_text", side_effect=AssertionError("Unexpected config read")):
                main(self.options + ["--check-only", "--first-time", "10"])
        self.assertIn("Layout: (4, 3, 2, 3)", output.getvalue())
        self.assertIn("first/last time: 10/10.06 s", output.getvalue())
        self.assertFalse(Path("Animations_velocity").exists())

    def test_cli_requires_explicit_file_dimensions_and_dt(self):
        with patch.dict(os.environ, {}, clear=True), contextlib.redirect_stderr(io.StringIO()):
            for name in ("--plane-file", "--nx", "--ny", "--dt"):
                options = self.options.copy()
                index = options.index(name)
                del options[index:index + 2]
                with self.assertRaises(SystemExit) as error:
                    main(options + ["--check-only"])
                self.assertEqual(error.exception.code, 2)
            options = self.options.copy()
            options[options.index("--nx") + 1] = "3.5"
            with self.assertRaises(SystemExit) as error:
                main(options + ["--check-only"])
            self.assertEqual(error.exception.code, 2)

    def run_renderer_process(self, code):
        """Include interpreter shutdown so PyGMT's atexit errors are observable."""
        environment = os.environ.copy()
        for key in ("OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE", "PMI_RANK", "PMI_SIZE",
                    "PMIX_RANK", "PMIX_SIZE", "SLURM_PROCID", "SLURM_NTASKS"):
            environment.pop(key, None)
        environment["PYTHONPATH"] = str(Path(__file__).resolve().parents[1])
        scratch = Path("scratch").resolve()
        scratch.mkdir(exist_ok=True)
        environment["SLURM_TMPDIR"] = str(scratch)
        result = subprocess.run([sys.executable, "-c", code], env=environment,
                                capture_output=True, text=True, timeout=180)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertNotRegex(result.stderr, r"\[WARNING\]|\[ERROR\]|Exception ignored|Traceback")
        self.assertEqual(list(scratch.iterdir()), [])

    @unittest.skipUnless(os.environ.get("HERCULES_TEST_PYGMT") == "1",
                         "Set HERCULES_TEST_PYGMT=1 with PyGMT/GMT, xarray and joblib installed")
    def test_render_and_shutdown_with_serial_and_reused_parallel_workers(self):
        for workers in (1, 2):
            with self.subTest(workers=workers):
                output = Path(f"pngs_{workers}")
                options = self.options + ["--workers", str(workers), "--dpi", "50",
                                           "--output-dir", str(output)]
                self.run_renderer_process(f"""
import os
import tempfile
from pathlib import Path
from script.animations_mpi import main
tempfile.gettempdir()  # Joblib also initializes this cache in the parent.
def state():
    return (Path.cwd(), tempfile.tempdir,
            [os.environ.get(key) for key in ('GMT_USERDIR', 'GMT_TMPDIR', 'TMPDIR')])
before = state()
main({options!r})
main({options!r})  # Exercise another chunk in a process that has imported PyGMT.
assert state() == before
""")
                images = sorted(output.glob("*.png"))
                self.assertEqual(len(images), 4)
                for image in images:
                    self.assertTrue(image.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"))

    @unittest.skipUnless(os.environ.get("HERCULES_TEST_PYGMT") == "1",
                         "Set HERCULES_TEST_PYGMT=1 with PyGMT/GMT, xarray and joblib installed")
    def test_failed_frame_closes_session_and_allows_next_render(self):
        bad = np.fromfile(self.plane_file, dtype=np.float64).reshape(4, 3, 2, 3)
        bad[1, 0, 0, 2] = np.nan
        bad.tofile("planes/bad.0")
        good_options = self.options + ["--workers", "1", "--dpi", "50"]
        bad_options = good_options.copy()
        bad_options[bad_options.index("--plane-file") + 1] = "planes/bad.0"
        self.run_renderer_process(f"""
import os
from pathlib import Path
from script.animations_mpi import main
before = (Path.cwd(), dict(os.environ))
try:
    main({bad_options!r})
except ValueError as error:
    assert 'Non-finite velocity' in str(error)
else:
    raise AssertionError('Invalid frame must fail')
assert (Path.cwd(), dict(os.environ)) == before
assert not list(Path(os.environ['SLURM_TMPDIR']).iterdir())
main({good_options!r})
""")


if __name__ == "__main__":
    unittest.main()
