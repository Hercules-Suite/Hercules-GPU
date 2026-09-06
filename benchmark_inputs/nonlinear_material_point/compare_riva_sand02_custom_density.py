"""Compare native custom-bound diagnostics to an untouched private oracle.

Usage: python3 compare_riva_sand02_custom_density.py PORT_ROOT NATIVE_TEST
No expected output or golden files are created or modified.
"""
from __future__ import annotations

import csv
import importlib
import io
import math
from pathlib import Path
import subprocess
import sys


def main() -> None:
    port = Path(sys.argv[1]).resolve()
    executable = Path(sys.argv[2]).resolve()
    sys.path.insert(0, str(port.parent))
    model_module = importlib.import_module(f"{port.name}.model")
    reference = importlib.import_module(f"{port.name}.reference")
    import numpy as np

    result = subprocess.run([str(executable), "--probe"], check=True,
                            capture_output=True, text=True)
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    assert len(rows) == 15, f"expected 15 cases, got {len(rows)}"
    count = 0
    for row in rows:
        values = {key: float(value) for key, value in row.items()}
        model = model_module.RIVASandIntermediateBiasFlowModel(
            reference.reference_parameters(e_max=values["e_max"],
                                           e_min=values["e_min"]))
        stress = np.array([[-19.4, 10.0, 0.0], [10.0, -19.4, 0.0],
                           [0.0, 0.0, -40.0]])
        zero_shear = np.diag(np.diag(stress))
        model.initialize(stress, values["e0"])
        model.begin_cyclic_phase(reference_stress=zero_shear)
        state = model.state
        state.amplitude_reversals = 4
        state.cyclic_amplitude = 0.6
        state.cyclic_direction = np.array([[0.0, 1/math.sqrt(2.0), 0.0],
                                          [1/math.sqrt(2.0), 0.0, 0.0],
                                          [0.0, 0.0, 0.0]])
        pressure = -float(np.trace(state.stress))/3.0
        expected = {
            "reference_density": model.reference_relative_density(),
            "initial_density": model.initial_relative_density(state),
            "raw": model.raw_state_shakedown_factor(pressure, state),
            "shakedown": model.state_shakedown_factor(pressure, state),
            "shear": model.state_shakedown_shear_multiplier(pressure, state),
            "density_weight": model.dense_state_weight(state),
            "dense_bias_gate": model.phase_dense_bias_gate(state),
            "dense_transition": model.phase_volume_dense_transition(state),
            "branch_compliance": model.branch_compliance_multiplier(state),
            "phase_hardening": model.bias_hardening_boost(state),
        }
        for name, target in expected.items():
            assert math.isclose(values[name], target, rel_tol=1e-12, abs_tol=1e-13), (
                values["e_max"], values["e_min"], values["e0"],
                name, values[name], target)
            count += 1
    print(f"PASS: {len(rows)} custom-bound cases, {count} native/Python diagnostics")


if __name__ == "__main__":
    main()
