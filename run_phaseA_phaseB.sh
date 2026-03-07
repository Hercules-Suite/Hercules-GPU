#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CASE_DIR="${CASE_DIR:-$PWD}"

# Default runtime tuning knobs (can still be overridden by user env).
: "${OMP_PROC_BIND:=close}"
: "${OMP_PLACES:=cores}"
: "${HERCULES_ADJUST_ASSIGN_OMP:=1}"
: "${HERCULES_ADJUST_DIST_OMP:=1}"
: "${HERCULES_SPARSE_CLEAR_MODE:=auto}"
: "${HERCULES_STATION_IO_OMP:=1}"
: "${HERCULES_PLANE_IO_OMP:=1}"
export OMP_PROC_BIND OMP_PLACES
export HERCULES_ADJUST_ASSIGN_OMP HERCULES_ADJUST_DIST_OMP
export HERCULES_SPARSE_CLEAR_MODE HERCULES_STATION_IO_OMP HERCULES_PLANE_IO_OMP

MPI_RUN="${MPI_RUN:-mpirun}"
BASE_IN="${BASE_IN:-}"
if [[ -f "${CASE_DIR}/inputfiles/Test_A.in" ]]; then
  PHASE_A_DEFAULT_IN="${CASE_DIR}/inputfiles/Test_A.in"
else
  PHASE_A_DEFAULT_IN="${ROOT_DIR}/Input/Test_A.in"
fi
if [[ -f "${CASE_DIR}/inputfiles/Test_B.in" ]]; then
  PHASE_B_DEFAULT_IN="${CASE_DIR}/inputfiles/Test_B.in"
else
  PHASE_B_DEFAULT_IN="${ROOT_DIR}/Input/Test_B.in"
fi
PHASE_A_IN="${PHASE_A_IN:-${PHASE_A_DEFAULT_IN}}"
PHASE_B_IN="${PHASE_B_IN:-${PHASE_B_DEFAULT_IN}}"
PHASE_A_MPI_ARGS="${PHASE_A_MPI_ARGS:-}"
PHASE_B_MPI_ARGS="${PHASE_B_MPI_ARGS:-}"

if [[ $# -ge 3 ]]; then
  BASE_IN="$3"
fi

if [[ -n "${PSOLVE_BIN:-}" ]]; then
  PSOLVE_BIN="${PSOLVE_BIN}"
elif [[ -x "${SCRIPT_DIR}/quake/forward/psolve" ]]; then
  PSOLVE_BIN="${SCRIPT_DIR}/quake/forward/psolve"
elif [[ -x "${SCRIPT_DIR}/psolve" ]]; then
  PSOLVE_BIN="${SCRIPT_DIR}/psolve"
else
  echo "ERROR: could not find psolve binary. Set PSOLVE_BIN=/path/to/psolve" >&2
  exit 1
fi

TMP_PHASE_A_IN=""
TMP_PHASE_B_IN=""

make_phase_input() {
  local src="$1"
  local dst="$2"
  local out_mesh="$3"
  local mesh_only="$4"
  local mesh_input="$5"

  awk -v outm="${out_mesh}" -v monly="${mesh_only}" -v minput="${mesh_input}" '
    BEGIN {
      seen_out = 0;
      seen_only = 0;
      seen_input = 0;
    }
    /^[[:space:]]*output_mesh[[:space:]]*=/ {
      print "output_mesh                 = " outm;
      seen_out = 1;
      next;
    }
    /^[[:space:]]*mesh_only_mode[[:space:]]*=/ {
      print "mesh_only_mode              = " monly;
      seen_only = 1;
      next;
    }
    /^[[:space:]]*mesh_input_mode[[:space:]]*=/ {
      print "mesh_input_mode             = " minput;
      seen_input = 1;
      next;
    }
    { print }
    END {
      if (!seen_out) {
        print "output_mesh                 = " outm;
      }
      if (!seen_only) {
        print "mesh_only_mode              = " monly;
      }
      if (!seen_input) {
        print "mesh_input_mode             = " minput;
      }
    }
  ' "${src}" >"${dst}"
}

if [[ -n "${BASE_IN}" ]]; then
  if [[ ! -f "${BASE_IN}" ]]; then
    echo "ERROR: BASE_IN file not found: ${BASE_IN}" >&2
    exit 1
  fi

  TMPDIR_ROOT="${TMPDIR:-/tmp}"
  TMP_PHASE_A_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phaseA.XXXXXX.in")"
  TMP_PHASE_B_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phaseB.XXXXXX.in")"

  cleanup_tmp_inputs() {
    [[ -n "${TMP_PHASE_A_IN}" ]] && rm -f "${TMP_PHASE_A_IN}"
    [[ -n "${TMP_PHASE_B_IN}" ]] && rm -f "${TMP_PHASE_B_IN}"
  }
  trap cleanup_tmp_inputs EXIT

  make_phase_input "${BASE_IN}" "${TMP_PHASE_A_IN}" 1 1 0
  make_phase_input "${BASE_IN}" "${TMP_PHASE_B_IN}" 0 0 1
  PHASE_A_IN="${TMP_PHASE_A_IN}"
  PHASE_B_IN="${TMP_PHASE_B_IN}"
fi

detect_alloc_nodes() {
  if [[ -n "${SLURM_JOB_NUM_NODES:-}" ]]; then
    echo "${SLURM_JOB_NUM_NODES}"
    return
  fi
  if [[ -n "${PBS_NNODES:-}" ]]; then
    echo "${PBS_NNODES}"
    return
  fi
  if [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]]; then
    sort -u "${PBS_NODEFILE}" | wc -l | tr -d ' '
    return
  fi
  if [[ -n "${SLURM_NODELIST:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    scontrol show hostnames "${SLURM_NODELIST}" | wc -l | tr -d ' '
    return
  fi
  echo "0"
}

ALLOC_NODES="$(detect_alloc_nodes)"

if [[ $# -ge 1 ]]; then
  PHASE_A_NP="$1"
else
  PHASE_A_NP="${PHASE_A_NP:-}"
fi

if [[ $# -ge 2 ]]; then
  PHASE_B_NP="$2"
else
  if [[ "${ALLOC_NODES}" -gt 0 ]]; then
    PHASE_B_NP="${PHASE_B_NP:-${ALLOC_NODES}}"
  else
    PHASE_B_NP="${PHASE_B_NP:-1}"
  fi
fi

if [[ ! -f "${PHASE_A_IN}" ]]; then
  echo "ERROR: Phase A input file not found: ${PHASE_A_IN}" >&2
  exit 1
fi
if [[ ! -f "${PHASE_B_IN}" ]]; then
  echo "ERROR: Phase B input file not found: ${PHASE_B_IN}" >&2
  exit 1
fi

if [[ -n "${PHASE_A_NP}" ]] && ( ! [[ "${PHASE_A_NP}" =~ ^[0-9]+$ ]] || [[ "${PHASE_A_NP}" -lt 1 ]] ); then
  echo "ERROR: invalid Phase A rank count: ${PHASE_A_NP}" >&2
  exit 1
fi

if ! [[ "${PHASE_B_NP}" =~ ^[0-9]+$ ]] || [[ "${PHASE_B_NP}" -lt 1 ]]; then
  echo "ERROR: invalid Phase B rank count: ${PHASE_B_NP}" >&2
  exit 1
fi

if [[ "${ALLOC_NODES}" -gt 0 ]] && [[ "${PHASE_B_NP}" -gt "${ALLOC_NODES}" ]]; then
  echo "ERROR: requested Phase B ranks (${PHASE_B_NP}) exceed allocated nodes (${ALLOC_NODES})." >&2
  echo "       Use <= allocated node count to guarantee one rank per node." >&2
  exit 1
fi

echo "=== Two-phase run ==="
echo "MPI launcher:   ${MPI_RUN}"
echo "psolve binary:  ${PSOLVE_BIN}"
echo "Case dir:       ${CASE_DIR}"
if [[ -n "${BASE_IN}" ]]; then
  echo "Base input:     ${BASE_IN}"
  echo "Phase A input:  ${PHASE_A_IN} (generated)"
  echo "Phase B input:  ${PHASE_B_IN} (generated)"
else
  echo "Phase A input:  ${PHASE_A_IN}"
  echo "Phase B input:  ${PHASE_B_IN}"
fi
if [[ "${ALLOC_NODES}" -gt 0 ]]; then
  echo "Alloc nodes:    ${ALLOC_NODES}"
else
  echo "Alloc nodes:    unknown (scheduler env not detected)"
fi
if [[ -n "${PHASE_A_NP}" ]]; then
  echo "Phase A ranks:  ${PHASE_A_NP}"
else
  echo "Phase A ranks:  launcher default"
fi
echo "Phase B ranks:  ${PHASE_B_NP} (1 rank/node)"
echo

if [[ -n "${PHASE_A_NP}" ]]; then
  echo "[Phase A] ${MPI_RUN} -np ${PHASE_A_NP} ${PSOLVE_BIN} ${PHASE_A_IN}"
else
  echo "[Phase A] ${MPI_RUN} ${PSOLVE_BIN} ${PHASE_A_IN}"
fi
if [[ -n "${PHASE_A_MPI_ARGS}" ]]; then
  # shellcheck disable=SC2206
  A_ARGS=(${PHASE_A_MPI_ARGS})
else
  A_ARGS=()
fi
if [[ -n "${PHASE_A_NP}" ]]; then
  "${MPI_RUN}" -np "${PHASE_A_NP}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
else
  "${MPI_RUN}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
fi

echo
echo "[Phase B] ${MPI_RUN} -np ${PHASE_B_NP} --map-by ppr:1:node --rank-by node ${PSOLVE_BIN} ${PHASE_B_IN}"
if [[ -n "${PHASE_B_MPI_ARGS}" ]]; then
  # shellcheck disable=SC2206
  B_ARGS=(${PHASE_B_MPI_ARGS})
else
  B_ARGS=()
fi
"${MPI_RUN}" -np "${PHASE_B_NP}" --map-by ppr:1:node --rank-by node "${B_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_B_IN}"

echo
echo "Two-phase run completed."
