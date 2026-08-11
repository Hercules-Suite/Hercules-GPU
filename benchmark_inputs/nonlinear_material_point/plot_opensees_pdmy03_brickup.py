#!/usr/bin/env python3
"""Plot coupled OpenSees SSPbrickUP and Hercules PDMY03 BrickUP results."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from compare_opensees_pdmy03_brickup import (
    COMPONENTS,
    read_hercules,
    read_opensees,
)


TRIAXIAL_TEST = "pdmy03_sspbrickup_undrained_triaxial"
CYCLIC_TEST = "pdmy03_sspbrickup_undrained_cyclic_shear"


def effective_invariants(stress: list[float]) -> tuple[float, float]:
    """Return compression-positive p' and triaxial q in Pa."""
    mean = sum(stress[:3]) / 3.0
    deviator = [stress[index] - mean for index in range(3)]
    shear_square = sum(value * value for value in stress[3:])
    q = math.sqrt(
        1.5 * (sum(value * value for value in deviator) + 2.0 * shear_square)
    )
    return -mean, q


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot coupled Hercules PDMY03 BrickUP and OpenSees SSPbrickUP "
            "triaxial or cyclic-shear histories."
        )
    )
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("pdmy03_sspbrickup_hercules_opensees"),
    )
    args = parser.parse_args()

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required; activate a Python environment containing "
            "matplotlib and rerun this script"
        ) from exc

    hercules = read_hercules(args.hercules_csv)
    test_id, opensees = read_opensees(args.opensees_csv)
    if len(hercules) != len(opensees):
        raise SystemExit(
            "Hercules and OpenSees histories have different lengths: "
            f"{len(hercules)} and {len(opensees)}"
        )

    h_effective = [
        [row[f"sig_{component}"] for component in COMPONENTS]
        for row in hercules
    ]
    o_effective = [
        [1000.0 * float(row[f"sigma_{component}_kpa"])
         for component in COMPONENTS]
        for row in opensees
    ]
    h_total = [
        [row[f"total_sig_{component}"] for component in COMPONENTS]
        for row in hercules
    ]
    o_total = [
        [1000.0 * float(row[f"total_sigma_{component}_kpa"])
         for component in COMPONENTS]
        for row in opensees
    ]
    h_invariants = [effective_invariants(stress) for stress in h_effective]
    o_invariants = [effective_invariants(stress) for stress in o_effective]
    h_p_kpa = [value[0] / 1000.0 for value in h_invariants]
    o_p_kpa = [value[0] / 1000.0 for value in o_invariants]
    h_q_kpa = [value[1] / 1000.0 for value in h_invariants]
    o_q_kpa = [value[1] / 1000.0 for value in o_invariants]
    h_pressure_kpa = [row["pore_pressure_pa"] / 1000.0 for row in hercules]
    o_pressure_kpa = [float(row["pore_pressure_kpa"]) for row in opensees]
    target_kpa = [float(row["confinement_kpa"]) for row in opensees]
    steps = list(range(len(hercules)))

    if test_id == TRIAXIAL_TEST:
        initial_axial_strain = hercules[0]["eps_zz"]
        control = [
            -100.0 * (row["eps_zz"] - initial_axial_strain)
            for row in hercules
        ]
        h_response = h_q_kpa
        o_response = o_q_kpa
        response_xlabel = r"Axial compression, $-\Delta\varepsilon_{zz}$ (%)"
        response_ylabel = r"Deviatoric stress, $q$ (kPa)"
        response_title = "Deviatoric stress response"
        h_confinement = [
            -0.5 * (stress[0] + stress[1]) / 1000.0 for stress in h_total
        ]
        o_confinement = [
            -0.5 * (stress[0] + stress[1]) / 1000.0 for stress in o_total
        ]
        confinement_label = "radial"
        progress = control
        progress_xlabel = response_xlabel
        path_ylabel = r"Deviatoric stress, $q$ (kPa)"
        h_path_response = h_q_kpa
        o_path_response = o_q_kpa
        title = "PDMY03 sealed undrained triaxial: OpenSees vs Hercules"
    elif test_id == CYCLIC_TEST:
        initial_shear_strain = hercules[0]["eps_xz"]
        control = [
            200.0 * (row["eps_xz"] - initial_shear_strain)
            for row in hercules
        ]
        h_response = [stress[5] / 1000.0 for stress in h_effective]
        o_response = [stress[5] / 1000.0 for stress in o_effective]
        response_xlabel = r"Engineering shear strain, $\gamma_{xz}$ (%)"
        response_ylabel = r"Shear stress, $\tau_{xz}$ (kPa)"
        response_title = "Cyclic stress–strain response"
        h_confinement = [-stress[2] / 1000.0 for stress in h_total]
        o_confinement = [-stress[2] / 1000.0 for stress in o_total]
        confinement_label = "vertical"
        progress = steps
        progress_xlabel = "Coupled state"
        path_ylabel = r"Shear stress, $\tau_{xz}$ (kPa)"
        h_path_response = h_response
        o_path_response = o_response
        title = "PDMY03 sealed cyclic shear: OpenSees vs Hercules"
    else:
        raise SystemExit(f"unsupported coupled PDMY03 test: {test_id}")

    h_confinement_error = [
        value - target for value, target in zip(h_confinement, target_kpa)
    ]
    o_confinement_error = [
        value - target for value, target in zip(o_confinement, target_kpa)
    ]

    hercules_color = "#0072B2"
    opensees_color = "#D55E00"
    linewidth = 1.65
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    fig.suptitle(title, fontsize=15)

    ax = axes[0, 0]
    ax.plot(control, h_response, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(control, o_response, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(response_xlabel)
    ax.set_ylabel(response_ylabel)
    ax.set_title(response_title)
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(progress, h_pressure_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(progress, o_pressure_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(progress_xlabel)
    ax.set_ylabel(r"Pore pressure, $u$ (kPa)")
    ax.set_title("Pore-pressure response")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(h_p_kpa, h_path_response, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(o_p_kpa, o_path_response, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Effective mean stress, $p'$ (kPa)")
    ax.set_ylabel(path_ylabel)
    ax.set_title("Effective stress path")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.axhline(0.0, color="0.45", linewidth=0.8)
    ax.plot(progress, h_confinement_error, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(progress, o_confinement_error, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(progress_xlabel)
    ax.set_ylabel(f"Total {confinement_label} confinement error (kPa)")
    ax.set_title(f"Confinement drift (target {target_kpa[0]:g} kPa)")
    ax.legend(frameon=False)

    for ax in axes.flat:
        ax.grid(True, color="0.88", linewidth=0.7)
        ax.tick_params(direction="in", top=True, right=True)
        for spine in ax.spines.values():
            spine.set_linewidth(0.8)

    prefix = args.output_prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    png_path = prefix.with_suffix(".png")
    pdf_path = prefix.with_suffix(".pdf")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Coupled canonical test: {test_id}")
    print(f"PNG output: {png_path.resolve()}")
    print(f"PDF output: {pdf_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
