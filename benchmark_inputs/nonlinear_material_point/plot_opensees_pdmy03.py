#!/usr/bin/env python3
"""Plot the synthetic Hercules/OpenSees PDMY03 characterization."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


HERCULES_STRESS_FIELDS = (
    "sig_xx", "sig_yy", "sig_zz", "sig_xy", "sig_yz", "sig_xz",
)
OPENSEES_STRESS_FIELDS = (
    "sigma_xx", "sigma_yy", "sigma_zz", "sigma_xy", "sigma_yz", "sigma_xz",
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise SystemExit(f"CSV has no header: {path}")
        rows = list(reader)
    if not rows:
        raise SystemExit(f"CSV is empty: {path}")
    return rows


def require_fields(
    path: Path, rows: list[dict[str, str]], fields: tuple[str, ...]
) -> None:
    missing = [field for field in fields if field not in rows[0]]
    if missing:
        raise SystemExit(f"{path} is missing fields: {', '.join(missing)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot Hercules and OpenSees PDMY03 cyclic-shear results."
    )
    parser.add_argument("--hercules-csv", type=Path, required=True)
    parser.add_argument("--opensees-csv", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("opensees_pdmy03_comparison"),
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

    hercules = read_rows(args.hercules_csv)
    opensees = read_rows(args.opensees_csv)
    require_fields(
        args.hercules_csv,
        hercules,
        ("eps_xy",) + HERCULES_STRESS_FIELDS,
    )
    require_fields(
        args.opensees_csv,
        opensees,
        ("test_id", "driver_version", "gamma_xy") + OPENSEES_STRESS_FIELDS,
    )
    if len(hercules) != len(opensees):
        raise SystemExit(
            "Hercules and OpenSees histories have different lengths: "
            f"{len(hercules)} and {len(opensees)}"
        )
    if any(
        row["test_id"] != "pdmy03_cyclic_shear"
        or row["driver_version"] != "1"
        for row in opensees
    ):
        raise SystemExit(
            "OpenSees CSV was not produced by the PDMY03 cyclic-shear "
            "driver version 1"
        )

    steps = list(range(len(hercules)))
    gamma_h_percent = [200.0 * float(row["eps_xy"]) for row in hercules]
    gamma_o_percent = [100.0 * float(row["gamma_xy"]) for row in opensees]
    stress_h = {
        field: [float(row[field]) for row in hercules]
        for field in HERCULES_STRESS_FIELDS
    }
    stress_o = {
        field: [float(row[field]) for row in opensees]
        for field in OPENSEES_STRESS_FIELDS
    }
    tau_h_kpa = [value / 1000.0 for value in stress_h["sig_xy"]]
    tau_o_kpa = [value / 1000.0 for value in stress_o["sigma_xy"]]
    p_h_kpa = [
        -(sxx + syy + szz) / 3000.0
        for sxx, syy, szz in zip(
            stress_h["sig_xx"], stress_h["sig_yy"], stress_h["sig_zz"]
        )
    ]
    p_o_kpa = [
        -(sxx + syy + szz) / 3000.0
        for sxx, syy, szz in zip(
            stress_o["sigma_xx"], stress_o["sigma_yy"], stress_o["sigma_zz"]
        )
    ]
    maximum_stress_difference = []
    shear_stress_difference = []
    for index in steps:
        differences = [
            abs(stress_h[hfield][index] - stress_o[ofield][index])
            for hfield, ofield in zip(
                HERCULES_STRESS_FIELDS, OPENSEES_STRESS_FIELDS
            )
        ]
        maximum_stress_difference.append(max(max(differences), 1.0e-16))
        shear_stress_difference.append(max(differences[3], 1.0e-16))

    hercules_color = "#0072B2"
    opensees_color = "#D55E00"
    difference_color = "#009E73"
    linewidth = 1.7

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    fig.suptitle("PressureDependMultiYield03: Hercules vs OpenSees", fontsize=15)

    ax = axes[0, 0]
    ax.plot(gamma_h_percent, tau_h_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(gamma_o_percent, tau_o_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Engineering shear strain, $\gamma_{xy}$ (%)")
    ax.set_ylabel(r"Shear stress, $\tau_{xy}$ (kPa)")
    ax.set_title("Cyclic stress–strain response")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(gamma_h_percent, p_h_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(gamma_o_percent, p_o_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Engineering shear strain, $\gamma_{xy}$ (%)")
    ax.set_ylabel(r"Effective mean stress, $p'$ (kPa)")
    ax.set_title("Pressure-dependent response")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(steps, gamma_h_percent, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(steps, gamma_o_percent, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel("Material-point state")
    ax.set_ylabel(r"Engineering shear strain, $\gamma_{xy}$ (%)")
    ax.set_title("Prescribed strain history")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.semilogy(steps, maximum_stress_difference, color=difference_color,
                linewidth=linewidth, label=r"max $|\Delta\sigma_{ij}|$")
    ax.semilogy(steps, shear_stress_difference, color=opensees_color,
                linewidth=linewidth, linestyle="--",
                label=r"$|\Delta\tau_{xy}|$")
    ax.set_xlabel("Material-point state")
    ax.set_ylabel("Absolute stress difference (Pa)")
    ax.set_title("Cross-code difference")
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
