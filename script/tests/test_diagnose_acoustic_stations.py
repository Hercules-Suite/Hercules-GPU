import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "diagnose_acoustic_stations.py"
spec = importlib.util.spec_from_file_location("acoustic_diagnostics", SCRIPT)
diagnostics = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostics)


class AcousticStationDiagnosticsTest(unittest.TestCase):
    def test_constant_acceleration_and_positive_down_surface_sign(self):
        time = np.linspace(0., 4., 401)
        accel = np.array([.2, -.4, .1])
        u = .5 * time[:, None]**2 * accel
        v = time[:, None] * accel
        fields = dict(u=u, v=v, a=np.broadcast_to(accel, u.shape), p=-1025*9.81*u[:, 2:3])
        result = diagnostics.metrics(time, fields, 3., 1025., True)
        np.testing.assert_allclose(result["u_integral_v_max_residual_m"], 0., atol=1e-12)
        np.testing.assert_allclose(result["v_integral_a_max_residual_m_s"], 0., atol=1e-12)
        self.assertLess(result["surface_relation_rms_m"], 1e-12)
        fields["u"] = u + np.column_stack((time * .1, time * 0, time * 0))
        drift = diagnostics.metrics(time, fields, 3., 1025.)
        self.assertAlmostEqual(drift["u_integral_v_max_residual_m"][0], .4)

    def test_header_controls_column_order_and_checks_completion(self):
        names = ["Time(s)", "Pressure(Pa)"] + [c for k in ("v", "u", "a") for c in diagnostics.FIELDS[k]]
        data = np.arange(55., dtype=float).reshape(5, 11)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "station.18"
            np.savetxt(path, data, header=" ".join(names))
            time, fields = diagnostics.read_station(path, 5)
            np.testing.assert_array_equal(time, data[:, 0])
            np.testing.assert_array_equal(fields["u"], data[:, 5:8])
            np.testing.assert_array_equal(fields["p"], data[:, 1:2])
            with self.assertRaisesRegex(ValueError, "expected 6"):
                diagnostics.read_station(path, 6)
            data[2, 1] = np.nan
            np.savetxt(path, data, header=" ".join(names))
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                diagnostics.read_station(path)

    def test_owner_and_position_must_match(self):
        probe = dict(station=18, x_m=500, y_m=23250, z_m=2812.5)
        medium = dict(owner_medium="elastic", pressure_defined="no", rho_kgm3="1025",
                      x_m="500", y_m="23250", z_m="2812.5")
        with self.assertRaisesRegex(ValueError, "acoustic owner"):
            diagnostics.check_probe(probe, medium)
        medium.update(owner_medium="acoustic", pressure_defined="yes", z_m="0")
        with self.assertRaisesRegex(ValueError, "z_m disagrees"):
            diagnostics.check_probe(probe, medium)
        medium["z_m"] = "2812.5"
        self.assertEqual(diagnostics.check_probe(probe, medium), 1025)


if __name__ == "__main__":
    unittest.main()
