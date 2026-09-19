"""Analytical checks for plane layout, velocity, and MPI frame ownership."""

import contextlib
import io
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from animations_mpi import (detect_rank_info, load_plane, main, partition,
                            scalar, vertical_velocity)


class PlaneAnimationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.previous_cwd = Path.cwd()
        os.chdir(self.temporary.name)
        Path("planes").mkdir()
        self.parameters = Path("parameters.in")
        self.parameters.write_text(
            "number_output_planes = 1\n"
            "simulation_delta_time_sec = 0.002\n"
            "output_planes_print_rate = 10\n"
            "output_planes_directory = planes\n"
            "output_planes =\n"
            "# lat lon depth ds ns dd nd strike dip\n"
            "0 0 0 500.0 3 325.0 2 0 0\n"
        )
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

    def test_binary_point_order_and_plot_coordinates(self):
        plane = load_plane(self.parameters)
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
        plane = load_plane(self.parameters)
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

    def test_separate_plane_file_and_selected_plane(self):
        Path("geometry.in").write_text("output_planes =\n"
                                       "0 0 0 100 2 100 2 90 0\n"
                                       "0 0 0 500 3 325 2 360 0\n")
        text = self.parameters.read_text().replace("number_output_planes = 1",
                                                   "number_output_planes = 2")
        self.parameters.write_text(text + "planes_input_file = geometry.in\n")
        plane = load_plane(self.parameters, 1, "planes/planedisplacements.0")
        self.assertEqual((plane.nx, plane.ny, plane.dt), (3, 2, 0.02))
        with self.assertRaisesRegex(ValueError, "requires strike=0"):
            load_plane(self.parameters, 0, "planes/planedisplacements.0")

    def test_incomplete_bytes_and_insufficient_snapshots(self):
        path = Path("planes/planedisplacements.0")
        original = path.read_bytes()
        for extra in (b"x", b"12345678"):
            path.write_bytes(original + extra)
            with self.assertRaisesRegex(ValueError, "not a multiple"):
                load_plane(self.parameters)
        for size in (0, 3 * 2 * 3 * 8):
            path.write_bytes(original[:size])
            with self.assertRaisesRegex(ValueError, "At least two"):
                load_plane(self.parameters)

    def test_invalid_sampling_and_nonfinite_velocity(self):
        original = self.parameters.read_text()
        for value in ("0", "-0.001", "nan", "inf"):
            self.parameters.write_text(original.replace("= 0.002", f"= {value}"))
            with self.assertRaisesRegex(ValueError, "finite and positive"):
                load_plane(self.parameters)
        self.parameters.write_text(original.replace("rate = 10", "rate = 0"))
        with self.assertRaisesRegex(ValueError, "positive"):
            load_plane(self.parameters)
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

    def test_check_only_needs_no_plotting_libraries_and_restart_needs_time(self):
        output = io.StringIO()
        with patch.dict(os.environ, {}, clear=True), contextlib.redirect_stdout(output):
            main(["--check-only"])
            self.parameters.write_text(self.parameters.read_text() + "use_checkpoint = 1\n")
            with self.assertRaisesRegex(ValueError, "--first-time"):
                main(["--check-only"])
            main(["--check-only", "--first-time", "10"])
        self.assertIn("Layout: (4, 3, 2, 3)", output.getvalue())
        self.assertIn("first/last time: 10/10.06 s", output.getvalue())
        self.assertFalse(Path("Animations_velocity").exists())

    def test_parameter_names_are_matched_exactly(self):
        self.assertEqual(scalar(["output_planes_directory = planes",
                                 "output_planes_print_rate = 10"], "output_planes_print_rate"), "10")
        with self.assertRaisesRegex(ValueError, "Missing required"):
            scalar(["output_planes_print_rate_extra = 10"], "output_planes_print_rate")


if __name__ == "__main__":
    unittest.main()
