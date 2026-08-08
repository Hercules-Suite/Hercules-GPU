#!/usr/bin/env python3
"""Plot the canonical undrained Manzari--Dafalias comparison."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_FIELDS = (
    "step",
    "eps_xy",
    "hercules_sig_xy_pa",
    "opensees_sig_xy_pa",
    "hercules_p_eff_pa",
    "opensees_p_eff_pa",
    "hercules_q_pa",
    "opensees_q_pa",
    "hercules_ru",
    "opensees_ru",
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
    if not data["step"]:
        raise SystemExit(f"comparison CSV is empty: {path}")
    return data


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot Hercules and OpenSees stress, pore-pressure ratio, and "
            "effective-stress paths for the canonical undrained cyclic test."
        )
    )
    parser.add_argument("--comparison-csv", type=Path, required=True)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("manzari_md_undrained_hercules_opensees"),
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
    step = data["step"]
    gamma_percent = [200.0 * value for value in data["eps_xy"]]
    hercules_tau_kpa = [value / 1000.0
                        for value in data["hercules_sig_xy_pa"]]
    opensees_tau_kpa = [value / 1000.0
                        for value in data["opensees_sig_xy_pa"]]
    hercules_p_kpa = [value / 1000.0
                      for value in data["hercules_p_eff_pa"]]
    opensees_p_kpa = [value / 1000.0
                      for value in data["opensees_p_eff_pa"]]
    hercules_q_kpa = [value / 1000.0 for value in data["hercules_q_pa"]]
    opensees_q_kpa = [value / 1000.0 for value in data["opensees_q_pa"]]

    hercules_color = "#0072B2"
    opensees_color = "#D55E00"
    linewidth = 1.65

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    fig.suptitle("Undrained cyclic Manzari–Dafalias comparison", fontsize=15)

    ax = axes[0, 0]
    ax.plot(gamma_percent, hercules_tau_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(gamma_percent, opensees_tau_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Engineering shear strain, $\gamma_{xy}$ (%)")
    ax.set_ylabel(r"Shear stress, $\tau_{xy}$ (kPa)")
    ax.set_title("Cyclic stress–strain response")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    ax.plot(step, data["hercules_ru"], color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(step, data["opensees_ru"], color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel("Cyclic history state")
    ax.set_ylabel(r"Pore-pressure ratio, $r_u$")
    ax.set_title("Inferred pore-pressure buildup")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(step, hercules_p_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(step, opensees_p_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel("Cyclic history state")
    ax.set_ylabel(r"Effective mean stress, $p'$ (kPa)")
    ax.set_title("Effective-stress reduction")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    ax.plot(hercules_p_kpa, hercules_q_kpa, color=hercules_color,
            linewidth=linewidth, label="Hercules")
    ax.plot(opensees_p_kpa, opensees_q_kpa, color=opensees_color,
            linewidth=linewidth, linestyle="--", label="OpenSees")
    ax.set_xlabel(r"Effective mean stress, $p'$ (kPa)")
    ax.set_ylabel(r"Deviatoric stress, $q$ (kPa)")
    ax.set_title("Effective stress path")
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
