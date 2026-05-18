#!/bin/sh
# CPS analytical run for the ATTENUATING homogeneous halfspace
#   Vp=2000 Vs=1000 rho=2000 Qp=200 Qs=100  (constant Q)
# matching the Hercules type_of_damping=bkt2 model line
#   0.0  2000.0 1000.0 2000.0 200.0 100.0
# Same 5 offsets (0,4,8,12,16 km E of source) and 5 receiver depths,
# same 2-s quadratic-rate source as the no-atten reference -- the ONLY
# difference vs run_elastic_halfspace_no_atten_dipslip_5offsets.sh is Q,
# so v4 (bkt2 + cpml_pml_material_damping=bkt) can be validated against a
# matching-attenuation analytic instead of the no-atten sanity check.
#
# NOTE: CPS uses a causal frequency-independent (ETA=0) Q operator while
# Hercules BKT2 is a few-SLS band-fit to constant Q -- expect high
# correlation with a small dispersion-operator offset, not a bit match
# (same situation as the ~0.99 elastic comparison).
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

MODEL="elastic_halfspace_atten.mod"
DFILE="dfile_elastic_pml_5offsets"
SOURCE_DEPTH_KM="5.0"
RECEIVER_DEPTHS_KM="0.0000 0.5000 0.9175 0.9575 1.5000"

STRIKE_DEG="0"
DIP_DEG="90"
RAKE_DEG="90"
M0_DYNE_CM="1.0e22"
AZIMUTH_DEG="90"
SOURCE_FILE="../../../../hercules_discrete_source/cps_quadratic_rate_2s.txt"

OUT_BASE="runs/elastic_halfspace_atten_dipslip_5offsets/displacement"
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

echo "Attenuating halfspace 5-offsets CPS run is in ${OUT_BASE}"
