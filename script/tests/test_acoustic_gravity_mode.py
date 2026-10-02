"""Independent gravity dispersion and surface-gradient convergence checks."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from check_acoustic_gravity_mode import analytical_omega2, mode, production_element


class GravityModeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.element = production_element(Path(__file__).resolve().parents[2])

    def test_production_operator_preserves_constant_and_positive_energy(self):
        stiffness, _, _ = self.element
        np.testing.assert_allclose(stiffness.sum(axis=1), 0, atol=1e-15)
        eigenvalues = np.linalg.eigvalsh(stiffness)
        self.assertLess(abs(eigenvalues[0]), 1e-14)
        self.assertTrue(np.all(eigenvalues[1:] > 0))

    def test_gravity_frequency_converges_quadratically(self):
        rows = [mode(*self.element, 500., 175.78125, n) for n in (8, 16, 32)]
        errors = [abs(r['frequency_relative_error']) for r in rows]
        self.assertLess(errors[-1], 2e-4)
        for coarse, fine in zip(errors, errors[1:]):
            self.assertGreater(coarse/fine, 3.8)
            self.assertLess(coarse/fine, 4.2)
        self.assertTrue(all(r['eigen_residual_relative'] < 1e-7 for r in rows))

    def test_recovered_surface_motion_converges_more_slowly(self):
        rows = [mode(*self.element, 500., 175.78125, n) for n in (8, 16, 32, 64)]
        errors = [abs(r['surface_gradient_motion_ratio']-1) for r in rows]
        self.assertTrue(all(a > b for a, b in zip(errors, errors[1:])))
        self.assertLess(errors[-1], .02)
        self.assertGreater(errors[0], .10)

    def test_compressible_relation_has_incompressible_limit(self):
        k, depth = 2*np.pi/1000, 175.78125
        exact = 9.81*k*np.tanh(k*depth)
        self.assertAlmostEqual(analytical_omega2(k, depth, 1e9)/exact, 1., places=12)
        self.assertLess(analytical_omega2(k, depth, 1500), exact)

    def test_permanent_bottom_motion_exposes_surface_recovery_offset(self):
        coarse = mode(*self.element, 2000., 175.78125, 8)
        fine = mode(*self.element, 2000., 175.78125, 32)
        key = 'unit_bottom_static_surface_recovered_uz_m'
        self.assertGreater(coarse[key], .05)
        self.assertLess(fine[key], coarse[key]/3.8)
        self.assertAlmostEqual(
            fine['unit_bottom_static_surface_recovered_ux_amplitude_m'] /
            fine['unit_bottom_static_surface_exact_ux_amplitude_m'], 1., places=3)


if __name__ == '__main__':
    unittest.main()
