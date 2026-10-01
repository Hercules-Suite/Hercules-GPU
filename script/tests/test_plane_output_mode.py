"""Run the production plane backend selector with a host C++ compiler.

No CUDA/MPI installation is needed: MPI_Abort is stubbed. This checks setup
selection and diagnostics only, not interpolation or distributed file writes.
Run from the repository root: python3 -m unittest discover -s script/tests -p test_plane_output_mode.py
"""

import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest


class PlaneOutputModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = (Path(__file__).resolve().parents[2] /
                  "quake/forward/io_planes.cc").read_text()
        start = source.index(
            "static void planes_select_output_mode(int32_t myID, int io_pool_pe_count) {")
        end = source.index("\nstatic void ", start + 1)
        harness = r"""
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
static int comm_solver = 0;
static const int ERROR = -1;
static int MPI_Abort(int, int) { fputs("MPI_Abort called\n", stderr); return 0; }
static int g_planeGpuOutputEnabled, g_planeGpuRootEnabled, g_planeGpuMpiioEnabled;
static int g_planePrintSequence, g_planeMpiBatchCount, g_planeMpiBatchStartSequence;
static int64_t g_planeMpiBatchFlushCount;
static int g_planeHostOnlyPointCount;
""" + source[start:end] + r"""
int main(int argc, char** argv) {
    if (argc > 2) setenv("HERCULES_PLANE_IO_MODE", argv[2], 1);
    else unsetenv("HERCULES_PLANE_IO_MODE");
    planes_select_output_mode(0, atoi(argv[1]));
    printf("FLAGS %d %d %d\n", g_planeGpuOutputEnabled,
           g_planeGpuRootEnabled, g_planeGpuMpiioEnabled);
}
"""
        temp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(temp.cleanup)
        cpp = Path(temp.name) / "selector.cc"
        cls.binary = Path(temp.name) / "selector"
        cpp.write_text(harness)
        subprocess.run(shlex.split(os.environ.get("CXX", "c++")) +
                       ["-std=c++11", "-Wall", "-Wextra", "-Werror", str(cpp),
                        "-o", str(cls.binary)], check=True)

    def run_mode(self, mode, pool=0):
        args = [str(self.binary), str(pool)]
        if mode is not None:
            args.append(mode)
        return subprocess.run(args, capture_output=True, text=True, check=False)

    def test_unset_empty_and_explicit_mpiio(self):
        for mode in (None, "", "mpiio", "MPIIO", "mpi-io"):
            with self.subTest(mode=mode):
                result = self.run_mode(mode)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("GPU interpolation + MPI-IO", result.stdout)
                self.assertIn("FLAGS 1 0 1", result.stdout)
                self.assertEqual(result.stderr, "")

    def test_explicit_alternatives_and_aliases(self):
        groups = [(("hybrid", "HYBRID", "gpu-root", "root", "gpu", "on", "yes", "1"),
                   "1 1 0", "GPU interpolation + rank-0 gather/fwrite"),
                  (("legacy", "LEGACY", "off", "no", "0"),
                   "0 0 0", "legacy CPU interpolation + rank-0/I/O-pool")]
        for modes, flags, diagnostic in groups:
            for mode in modes:
                with self.subTest(mode=mode):
                    result = self.run_mode(mode)
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertIn("FLAGS " + flags, result.stdout)
                    self.assertIn(diagnostic, result.stdout)
                    self.assertEqual(result.stderr, "")

    def test_dedicated_io_pool_preserves_legacy_fallback(self):
        for mode in (None, "", "mpiio", "hybrid", "legacy"):
            with self.subTest(mode=mode):
                result = self.run_mode(mode, pool=2)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("FLAGS 0 0 0", result.stdout)
                self.assertIn("legacy CPU interpolation + rank-0/I/O-pool", result.stdout)
                if mode == "legacy":
                    self.assertEqual(result.stderr, "")
                else:
                    self.assertIn("using the existing I/O-pool backend", result.stderr)

    def test_auto_and_unknown_values_abort_even_with_io_pool(self):
        for mode in ("auto", "AUTO", "invalid"):
            for pool in (0, 2):
                with self.subTest(mode=mode, pool=pool):
                    result = self.run_mode(mode, pool)
                    self.assertEqual(result.returncode, 1)
                    self.assertIn("MPI_Abort called", result.stderr)
                    self.assertIn("expected mpiio (default), hybrid, or legacy", result.stderr)
                    self.assertNotIn("Plane output backend:", result.stdout)


if __name__ == "__main__":
    unittest.main()
