#!/bin/sh
# CPS analytical reference for the 937.5 m WATER (acoustic) layer over an
# elastic halfspace, no attenuation -- used to validate the acoustic CPML.
#
# This mirrors run_elastic_halfspace_no_atten_dipslip_16stations.sh EXACTLY
# (the 16-station elastic test the user's Hercules run matches): strike 0,
# dip 90, rake 90, M0 1e22, azimuth 90, -ROT, source depth 5 km, dfile
# offsets 0/4/8/12 km. Only MODEL -> flat_water_elastic_no_atten.mod and
# the receiver depths (the four Hercules station depths) differ.
#
# *** SOURCE-TIME FUNCTION ***
# The Hercules flat-water input uses slipfunction_2s.in (2.0 s, 201 pts),
# CONFIRMED by the user, so the matching CPS STF is cps_quadratic_rate_2s.txt.
# (The other matched pair is the 0.5 s slipfunction.in <-> cps_quadratic_rate.txt
# used by the 16-station *elastic* script; do NOT use that here.)
# A mismatched STF silently invalidates the comparison.
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
# Hercules station depths: 0 m (surface), 917.5 m (IN the 937.5 m water
# column), 957.5 m and 1500 m (in the elastic halfspace).
RECEIVER_DEPTHS_KM="0.0000 0.9175 0.9575 1.5000"

STRIKE_DEG="0"
DIP_DEG="90"
RAKE_DEG="90"
M0_DYNE_CM="1.0e22"
AZIMUTH_DEG="90"
# 2.0 s STF -- matches the Hercules slipfunction_2s.in (user-confirmed).
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
