#!/bin/sh
# CPS analytical reference for the 937.5 m WATER (acoustic) layer over an
# elastic halfspace, no attenuation -- used to validate the acoustic CPML.
#
# IMPORTANT: this MUST use the SAME source as the Hercules run and as the
# elastic CPS reference (run_elastic_halfspace_no_atten_dipslip_5offsets.sh),
# otherwise the comparison is meaningless. The previously generated water
# reference (runs/hercules_match_no_atten, from run_moment_source_no_atten.sh)
# used rake=0 + a triangular hpulse96 -t -l 25 pulse and is NOT comparable to
# the Hercules dip-slip / 2 s run. This script mirrors the elastic 5-offsets
# generator exactly:
#   strike 0, dip 90, RAKE 90, M0 1e22, azimuth 90, -ROT,
#   STF = the 2 s quadratic slip-rate (cps_quadratic_rate_2s.txt) with -STEP,
#   source depth 5 km.
# Only MODEL, the receiver depths, the dfile, and the output dir differ.
set -eu

need() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing CPS command: $1" >&2
    exit 127
  fi
}

need hprep96
need hspec96
need hpulse96
need fmech96
need f96tosac

MODEL="flat_water_elastic_no_atten.mod"
DFILE="dfile_flat_water"
SOURCE_DEPTH_KM="5.0"
# Hercules station depths: 0 m (surface), 917.5 m (IN the water column),
# 957.5 m (just below the 937.5 m water/solid interface, in the elastic
# halfspace), 1500 m (deeper in the elastic halfspace).
RECEIVER_DEPTHS_KM="0.0000 0.9175 0.9575 1.5000"

STRIKE_DEG="0"
DIP_DEG="90"
RAKE_DEG="90"
M0_DYNE_CM="1.0e22"
AZIMUTH_DEG="90"
SOURCE_FILE="../../../../hercules_discrete_source/cps_quadratic_rate_2s.txt"

OUT_BASE="runs/flat_water_no_atten_dipslip/displacement"
mkdir -p "${OUT_BASE}"

for HR in ${RECEIVER_DEPTHS_KM}; do
  RUN_DIR="${OUT_BASE}/hr_${HR}"
  mkdir -p "${RUN_DIR}"
  cp "${MODEL}" "${RUN_DIR}/"
  cp "${DFILE}" "${RUN_DIR}/"

  (
    cd "${RUN_DIR}"
    hprep96 -M "${MODEL}" -d "${DFILE}" -HS "${SOURCE_DEPTH_KM}" -HR "${HR}"
    hspec96
    hpulse96 -D -F "${SOURCE_FILE}" -STEP | \
      fmech96 -S "${STRIKE_DEG}" -D "${DIP_DEG}" -R "${RAKE_DEG}" \
        -M0 "${M0_DYNE_CM}" -A "${AZIMUTH_DEG}" -ROT | \
      f96tosac -FMT 2
  )
done

echo "Flat-water (acoustic) no-attenuation dip-slip CPS run is in ${OUT_BASE}"
