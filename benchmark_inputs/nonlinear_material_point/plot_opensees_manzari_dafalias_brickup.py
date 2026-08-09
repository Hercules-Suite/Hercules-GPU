#!/usr/bin/env python3
"""Plot the coupled SSPbrickUP/Hercules BrickUP comparison."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_FIELDS = (
    "eps_zz",
    "hercules_p_eff_pa",
    "opensees_p_eff_pa",
    "hercules_q_pa",
    "opensees_q_pa",
    "hercules_pore_pressure_pa",
    "opensees_pore_pressure_pa",
    "hercules_total_sig_xx_pa",
    "hercules_total_sig_yy_pa",
    "opensees_total_sig_xx_pa",
    "opensees_total_sig_yy_pa",
)


def read_comparison(path: Path) -> dict[str, list[float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise SystemExit(f"comparison CSV has no header: {path}")
        missing = [field for field in REQUIRED_FIELDS
                   if field not in reader.fieldnames]
        if missing:
            raise SystemExit("comparison CSV is missing fields: "
                             + ", ".join(missing))
        data = {field: [] for field in REQUIRED_FIELDS}
        for row in reader:
            for field in REQUIRED_FIELDS:
                data[field].append(float(row[field]))
    if not data["eps_zz"]:
        raise SystemExit(f"comparison CSV is empty: {path}")
    return data


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot Hercules BrickUP and OpenSees SSPbrickUP results for the "
            "sealed undrained triaxial benchmark."
        )
    )
    parser.add_argument("--comparison-csv", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("manzari_md_sspbrickup_hercules_opensees"),
    )
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit(
            "matplotlib is required; load a Python environment that provides "
            "matplotlib and rerun this script"
        ) from exc

    data = read_comparison(args.comparison_csv)
    initial_eps_zz = data["eps_zz"][0]
    axial_compression_percent = [
        -100.0 * (value - initial_eps_zz) for value in data["eps_zz"]
    ]
    hercules_p_kpa = [value / 1000.0
                      for value in data["hercules_p_eff_pa"]]
    opensees_p_kpa = [value / 1000.0
                      for value in data["opensees_p_eff_pa"]]
    hercules_q_kpa = [value / 1000.0
                      for value in data["hercules_q_pa"]]
    opensees_q_kpa = [value / 1000.0
                      for value in data["opensees_q_pa"]]
    hercules_pressure_kpa = [
        value / 1000.0 for value in data["hercules_pore_pressure_pa"]
    ]
    opensees_pressure_kpa = [
        value / 1000.0 for value in data["opensees_pore_pressure_pa"]
    ]
    hercules_radial_kpa = [
        -0.5 * (sxx + syy) / 1000.0
        for sxx, syy in zip(
            data["hercules_total_sig_xx_pa"],
            data["hercules_total_sig_yy_pa"],
        )
    ]
    opensees_radial_kpa = [
        -0.5 * (sxx + syy) / 1000.0
        for sxx, syy in zip(
            data["opensees_total_sig_xx_pa"],
            data["opensees_total_sig_yy_pa"],
        )
    ]
    target_radial_kpa = opensees_radial_kpa[0]
    hercules_radial_error_kpa = [
        value - target_radial_kpa for value in hercules_radial_kpa
    ]
    opensees_radial_error_kpa = [
        value - target_radial_kpa for value in opensees_radial_kpa
    ]

    hercules_color = "#0072B2"
    opensees_color = "#D55E00"
    linewidth = 1.65

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    fig.suptitle(
        "Sealed undrained triaxial: SSPbrickUP vs Hercules BrickUP",
        fontsize=15,
    )

    ax = axes[0, 0]
    ax.plot(axial_compression_percent, hercules_q_kpa,
            color=hercules_color, linewidth=linewidth, label="Hercules")
    ax.plot(axial_compression_percent, opensees_q_kpa,
            color=opensees_color, linewidth=linewidth, linestyle="--",
            label="OpenSees")
    ax.set_xlabel(r"Axial compression, $-\Delta\varepsilon_{zz}$ (%)")
    ax.set_ylabel(r"Deviatoric stress, $q$ (kPa)")
    ax.set_title("Deviatoric stress response")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(axial_compression_percent, hercules_pressure_kpa,
            color=hercules_color, linewidth=linewidth, label="Hercules")
    ax.plot(axial_compression_percent, opensees_pressure_kpa,
            color=opensees_color, linewidth=linewidth, linestyle="--",
            label="OpenSees")
    ax.set_xlabel(r"Axial compression, $-\Delta\varepsilon_{zz}$ (%)")
    ax.set_ylabel(r"Pore pressure, $u$ (kPa)")
    ax.set_title("Pore-pressure buildup")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(hercules_p_kpa, hercules_q_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(opensees_p_kpa, opensees_q_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Effective mean stress, $p'$ (kPa)")
    ax.set_ylabel(r"Deviatoric stress, $q$ (kPa)")
    ax.set_title("Effective stress path")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.axhline(0.0, color="0.45", linewidth=0.8)
    ax.plot(axial_compression_percent, hercules_radial_error_kpa,
            color=hercules_color, linewidth=linewidth, label="Hercules")
    ax.plot(axial_compression_percent, opensees_radial_error_kpa,
            color=opensees_color, linewidth=linewidth, linestyle="--",
            label="OpenSees")
    ax.set_xlabel(r"Axial compression, $-\Delta\varepsilon_{zz}$ (%)")
    ax.set_ylabel("Total radial confinement error (kPa)")
    ax.set_title(f"Confinement drift (target {target_radial_kpa:g} kPa)")
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

    print(f"PNG output: {png_path.resolve()}")
    print(f"PDF output: {pdf_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
