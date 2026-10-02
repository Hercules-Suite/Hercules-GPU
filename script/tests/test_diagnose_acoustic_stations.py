import importlib.util
import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

SCRIPT = Path(__file__).resolve().parents[1] / "diagnose_acoustic_stations.py"
spec = importlib.util.spec_from_file_location("acoustic_diagnostics", SCRIPT)
diagnostics = importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostics)


class AcousticStationDiagnosticsTest(unittest.TestCase):
    def test_changed_cpml_membership_excludes_surface_relation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            probes = root / 'probes.json'
            probes.write_text(json.dumps([dict(station=0, label='surface', surface=True,
                                               cpml=0, cpml_by_case={'thick': 1})]))
            names = ['Time(s)'] + [c for columns in diagnostics.FIELDS.values() for c in columns]
            time = np.linspace(0, 4, 41)
            data = np.zeros((len(time), len(names)))
            data[:, 0] = time
            for name in ('base', 'thick'):
                folder = root/name/'output/stations'
                folder.mkdir(parents=True)
                (folder/'station_mediums.txt').write_text(
                    '# station_id owner_medium pressure_defined rho_kgm3\n0 acoustic yes 1025\n')
                np.savetxt(folder/'station.0', data, header=' '.join(names))
            args = argparse.Namespace(cases=[f'base={root}/base', f'thick={root}/thick'],
                probes=probes, gravity_cases=['base', 'thick'], late_start=3,
                history_stride=None, output_dir=root/'report', expected_samples=len(time), plot=False)
            diagnostics.run(args)
            cases = json.loads((args.output_dir/'summary.json').read_text())['stations'][0]['cases']
            self.assertIn('surface_relation_rms_m', cases['base'])
            self.assertNotIn('surface_relation_rms_m', cases['thick'])
            self.assertEqual(cases['thick']['cpml'], 1)
            self.assertEqual(diagnostics.probe_cpml({'cpml': 0, 'cpml_by_case': {'thick': None}}, 'thick'), None)

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
        medium["z_m"] = "nan"
        with self.assertRaisesRegex(ValueError, "z_m disagrees"):
            diagnostics.check_probe(probe, medium)


if __name__ == "__main__":
    unittest.main()
