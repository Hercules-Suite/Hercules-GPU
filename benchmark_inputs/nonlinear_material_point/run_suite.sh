#!/usr/bin/env bash
set -uo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 /absolute/or/relative/path/to/psolve" >&2
  exit 2
fi

psolve_bin=$1
if [[ $psolve_bin != /* ]]; then
  if ! psolve_dir=$(cd "$(dirname "$psolve_bin")" && pwd); then
    echo "Cannot resolve psolve directory: $(dirname "$psolve_bin")" >&2
    exit 2
  fi
  psolve_bin="$psolve_dir/$(basename "$psolve_bin")"
fi
if [[ ! -x $psolve_bin ]]; then
  echo "psolve is not executable: $psolve_bin" >&2
  exit 2
fi

if ! suite_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); then
  echo "Cannot resolve material-point suite directory" >&2
  exit 2
fi
cd "$suite_dir"

fixtures=(
  vonmises_ep_synthetic.in
  vonmises_fa_quartic_synthetic.in
  vonmises_fam_quartic_reversal_synthetic.in
  vonmises_bae_virgin_onset_synthetic.in
  vonmises_bae_legacy_id2_synthetic.in
  vonmises_bah_synthetic.in
  vonmises_gqh_synthetic.in
  vonmises_mkz_synthetic.in
  vonmises_ro_synthetic.in
  manzari_dafalias_modified_euler_synthetic.in
  manzari_dafalias_forward_euler_synthetic.in
  manzari_dafalias_forward_euler_demand_synthetic.in
  manzari_dafalias_seeded_onset_synthetic.in
  manzari_dafalias_pressure_floor_synthetic.in
  pdmy03_cyclic_shear_synthetic.in
  pdmy03_zero_tension_onset_synthetic.in
)

validate_md_pressure_floor() {
  local csv=manzari_dafalias_pressure_floor_synthetic.cpu.csv
  awk -F, '
    function absolute(value) { return value < 0.0 ? -value : value }
    NR == 1 {
      for (column = 1; column <= NF; column++) {
        if ($column == "step") step_column = column
        if ($column == "sig_xx") xx_column = column
        if ($column == "sig_yy") yy_column = column
        if ($column == "sig_zz") zz_column = column
      }
      next
    }
    $step_column == 1 || $step_column == 2 {
      # CSV stress follows the Hercules tension-positive convention.
      pressure = -($xx_column + $yy_column + $zz_column)/3.0
      if (absolute(pressure - 10.13) > 1.0e-6) exit 1
      cutoff_states++
    }
    END {
      if (!step_column || !xx_column || !yy_column || !zz_column ||
          cutoff_states != 2) exit 1
    }
  ' "$csv"
}

validate_md_forward_euler_substeps() {
  local mode=$1
  local prefix=$2
  local cap=$3
  local backend
  for backend in cpu gpu; do
    local csv="${prefix}.${backend}.csv"
    if ! awk -F, -v mode="$mode" -v cap="$cap" '
      NR == 1 {
        for (column = 1; column <= NF; column++) {
          if ($column == "plastic_multiplier") multiplier_column = column
          if ($column == "accepted_substeps") substeps_column = column
        }
        next
      }
      $multiplier_column > 0.0 {
        plastic_steps++
        if ($substeps_column < 1 || $substeps_column > cap) exit 1
        if (mode == "fixed" && $substeps_column != cap) exit 1
        if (mode == "demand" && $substeps_column < cap) reductions++
      }
      END {
        if (!multiplier_column || !substeps_column || plastic_steps == 0)
          exit 1
        if (mode == "demand" && reductions == 0) exit 1
      }
    ' "$csv"; then
      echo "Forward-Euler ${mode} substep policy failed for ${csv}" >&2
      return 1
    fi
  done
}

failures=()
for fixture in "${fixtures[@]}"; do
  printf '\n==> %s\n' "$fixture"
  if ! mpirun -np 1 "$psolve_bin" --nonlinear-material-point "$fixture"; then
    failures+=("$fixture")
  elif [[ $fixture == manzari_dafalias_forward_euler_synthetic.in ]] &&
       ! validate_md_forward_euler_substeps fixed \
           manzari_dafalias_forward_euler_synthetic 100; then
    failures+=("$fixture")
  elif [[ $fixture == manzari_dafalias_forward_euler_demand_synthetic.in ]] &&
       ! validate_md_forward_euler_substeps demand \
           manzari_dafalias_forward_euler_demand_synthetic 100; then
    failures+=("$fixture")
  elif [[ $fixture == manzari_dafalias_pressure_floor_synthetic.in ]] &&
       ! validate_md_pressure_floor; then
    echo "Material-point pressure-floor states do not equal 1.0e-4 P_atm" >&2
    failures+=("$fixture")
  fi
done

if (( ${#failures[@]} > 0 )); then
  printf '\nFAIL: %d of %d nonlinear material-point fixtures failed:\n' \
    "${#failures[@]}" "${#fixtures[@]}" >&2
  printf '  %s\n' "${failures[@]}" >&2
  exit 1
fi

printf '\nPASS: all %d nonlinear material-point fixtures passed\n' \
  "${#fixtures[@]}"
