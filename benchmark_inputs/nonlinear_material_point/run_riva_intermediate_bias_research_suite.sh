#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 || $# > 2 )); then
  echo "Usage: $0 PORT_ROOT [CXX]" >&2
  exit 2
fi

port_root=$1
cxx=${2:-c++}
if [[ ! -d $port_root/golden_data/input_paths ||
      ! -f $port_root/verify_hercules_output.py ]]; then
  echo "PORT_ROOT is not a complete RIVA-Sand handoff package: $port_root" >&2
  exit 2
fi

if ! suite_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); then
  echo "Cannot resolve research suite directory" >&2
  exit 2
fi

suite_tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/riva-ib-suite.XXXXXX")
cleanup() {
  rm -rf -- "$suite_tmp_dir"
}
trap cleanup EXIT

compile_flags=(-std=c++17 -O2 -Wall -Wextra -pedantic)
driver_source=$suite_dir/riva_sand_intermediate_bias_research_native_driver.cc
state_source=$suite_dir/riva_sand_intermediate_bias_research_state_test.cc
benchmark_source=$suite_dir/riva_sand_intermediate_bias_research_benchmark.cc
driver_bin=$suite_tmp_dir/riva_ib_replay
state_bin=$suite_tmp_dir/riva_ib_state_test
benchmark_bin=$suite_tmp_dir/riva_ib_benchmark
native_output=$suite_tmp_dir/native_output
input_paths=$port_root/golden_data/input_paths
benchmark_repetitions=${RIVA_IB_BENCHMARK_REPETITIONS:-2000}

echo "==> Compile native replay driver with $cxx"
"$cxx" "${compile_flags[@]}" "$driver_source" -o "$driver_bin"

echo "==> Compile host state-contract test with $cxx"
"$cxx" "${compile_flags[@]}" "$state_source" -o "$state_bin"

echo "==> Compile research-vs-production benchmark with $cxx"
"$cxx" "${compile_flags[@]}" "$benchmark_source" -o "$benchmark_bin"

echo "==> Compile corrected branch-reversal/state/checkpoint checks"
"$cxx" "${compile_flags[@]}" "$suite_dir/riva_sand02_branch_reversal_test.cc" \
  -o "$suite_tmp_dir/riva_branch_test"
"$suite_tmp_dir/riva_branch_test" "$suite_tmp_dir/branch_results"

echo "==> Compile velocity-input and stiffness checks"
"$cxx" "${compile_flags[@]}" "$suite_dir/riva_sand_velocity_test.cc" \
  -o "$suite_tmp_dir/riva_sand_velocity_test"
"$suite_tmp_dir/riva_sand_velocity_test" \
  "$suite_dir/riva_sand_brickup_synthetic.in" \
  "$suite_dir/riva_sand_low_bias_synthetic.in" \
  "$suite_dir/riva_sand_ramped_geostatic_admission_synthetic.in" \
  "$suite_dir/riva_sand_custom_equivalence.in" \
  "$suite_dir/riva_sand_intermediate_bias_research_synthetic.in" \
  "$suite_dir/riva_sand_custom02_equivalence.in" \
  "$suite_dir/riva_sand_velocity_synthetic.in" \
  "$suite_dir/riva_sand02_velocity_synthetic.in"
echo "==> Compile RIVASand02 custom-interface state checks"
"$cxx" "${compile_flags[@]}" "$suite_dir/riva_sand_g0_test.cc" \
  -o "$suite_tmp_dir/riva_sand_g0_test"
"$suite_tmp_dir/riva_sand_g0_test"
"$cxx" "${compile_flags[@]}" "$suite_dir/riva_sand02_custom_test.cc" \
  -o "$suite_tmp_dir/riva_sand02_custom_test"
"$cxx" "${compile_flags[@]}" "$suite_dir/riva_sand02_custom_density_test.cc" \
  -o "$suite_tmp_dir/riva_sand02_custom_density_test"

echo "==> Replay frozen oracle with the explicit historical reference"
"$driver_bin" "$input_paths" "$native_output" --legacy-reversal-reference

echo "==> Replay the corrected default on all six input paths"
"$driver_bin" "$input_paths" "$suite_tmp_dir/corrected_output"

suite_status=0
if ! "$suite_tmp_dir/riva_sand02_custom_test"; then
  suite_status=1
fi
if ! "$suite_tmp_dir/riva_sand02_custom_density_test"; then
  suite_status=1
fi
if ! python3 "$suite_dir/compare_riva_sand02_custom_density.py" \
    "$port_root" "$suite_tmp_dir/riva_sand02_custom_density_test"; then
  suite_status=1
fi
echo "==> Verify historical-reference output: core (not corrected-mode acceptance)"
if ! python3 "$port_root/verify_hercules_output.py" \
    "$port_root/golden_data" "$native_output" --mode core; then
  suite_status=1
fi

echo "==> Verify historical-reference output: full (not corrected-mode acceptance)"
if ! python3 "$port_root/verify_hercules_output.py" \
    "$port_root/golden_data" "$native_output" --mode full; then
  suite_status=1
fi

echo "==> Run restart/state-isolation/BrickUP host contract"
if ! "$state_bin" "$input_paths"; then
  suite_status=1
fi

echo "==> Benchmark research against production ($benchmark_repetitions repetitions)"
if ! "$benchmark_bin" "$input_paths" "$benchmark_repetitions"; then
  suite_status=1
fi

if (( suite_status == 0 )); then
  echo "PASS: intermediate-bias research native suite"
else
  echo "FAIL: intermediate-bias research native suite" >&2
fi
exit "$suite_status"
