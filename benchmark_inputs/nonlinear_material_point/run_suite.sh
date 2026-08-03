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
  vonmises_bae_legacy_id2_synthetic.in
  vonmises_bah_synthetic.in
  vonmises_gqh_synthetic.in
  vonmises_mkz_synthetic.in
  vonmises_ro_synthetic.in
)

failures=()
for fixture in "${fixtures[@]}"; do
  printf '\n==> %s\n' "$fixture"
  if ! mpirun -np 1 "$psolve_bin" --nonlinear-material-point "$fixture"; then
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
