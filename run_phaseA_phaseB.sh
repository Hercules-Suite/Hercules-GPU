#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
RUN_DIR="${PWD}"

# Default runtime tuning knobs (can still be overridden by user env).
: "${OMP_PROC_BIND:=close}"
: "${OMP_PLACES:=cores}"
: "${HERCULES_ADJUST_ASSIGN_OMP:=1}"
: "${HERCULES_ADJUST_DIST_OMP:=1}"
: "${HERCULES_SPARSE_CLEAR_MODE:=auto}"
: "${HERCULES_SPARSE_HOST_OMP:=0}"
: "${HERCULES_STATION_IO_OMP:=1}"
: "${HERCULES_PLANE_IO_OMP:=1}"
export OMP_PROC_BIND OMP_PLACES
export HERCULES_ADJUST_ASSIGN_OMP HERCULES_ADJUST_DIST_OMP
export HERCULES_SPARSE_CLEAR_MODE HERCULES_SPARSE_HOST_OMP
export HERCULES_STATION_IO_OMP HERCULES_PLANE_IO_OMP

MPI_RUN="${MPI_RUN:-mpirun}"
BASE_IN="${BASE_IN:-}"
if [[ -f "${RUN_DIR}/inputfiles/Test_A.in" ]]; then
  PHASE_A_DEFAULT_IN="${RUN_DIR}/inputfiles/Test_A.in"
else
  PHASE_A_DEFAULT_IN="${ROOT_DIR}/Input/Test_A.in"
fi
if [[ -f "${RUN_DIR}/inputfiles/Test_B.in" ]]; then
  PHASE_B_DEFAULT_IN="${RUN_DIR}/inputfiles/Test_B.in"
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

trim_ws() {
  local s="$1"
  # Strip leading/trailing whitespace without external dependencies.
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "${s}"
}

extract_config_value() {
  local file="$1"
  local key="$2"
  awk -v k="${key}" '
    {
      line = $0
      sub(/#.*/, "", line)
      if (line ~ "^[[:space:]]*" k "[[:space:]]*=") {
        sub("^[[:space:]]*" k "[[:space:]]*=[[:space:]]*", "", line)
        print line
        exit
      }
    }
  ' "${file}"
}

set_or_add_param() {
  local src="$1"
  local dst="$2"
  local key="$3"
  local value="$4"

  awk -v k="${key}" -v v="${value}" '
    BEGIN { seen = 0 }
    {
      if ($0 ~ "^[[:space:]]*" k "[[:space:]]*=") {
        printf "%s = %s\n", k, v
        seen = 1
      } else {
        print $0
      }
    }
    END {
      if (!seen) {
        printf "%s = %s\n", k, v
      }
    }
  ' "${src}" >"${dst}"
}

resolve_mesh_path() {
  local p="$1"
  local run_dir="$2"
  local base_in_dir="$3"
  local run_path=""
  local base_path=""

  if [[ -z "${p}" ]]; then
    return
  fi

  if [[ "${p}" = /* ]]; then
    printf '%s' "${p}"
    return
  fi

  run_path="${run_dir}/${p}"
  base_path="${base_in_dir}/${p}"

  # Prefer run-directory resolution for generated two-phase runs.
  # If only an inputfile-relative file already exists, preserve compatibility.
  if [[ -e "${base_path}" && ! -e "${run_path}" ]]; then
    printf '%s' "${base_path}"
  else
    printf '%s' "${run_path}"
  fi
}

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

  # Normalize mesh etree path so both phases use the same absolute location.
  BASE_IN_DIR="$(cd -- "$(dirname -- "${BASE_IN}")" && pwd)"
  RUN_DIR_ABS="$(cd -- "${RUN_DIR}" && pwd)"
  MESH_ETREE_RAW="$(extract_config_value "${BASE_IN}" "mesh_etree_output_file")"
  MESH_ETREE_RAW="$(trim_ws "${MESH_ETREE_RAW}")"
  if [[ -n "${MESH_ETREE_RAW}" ]]; then
    EFFECTIVE_MESH_ETREE_PATH="$(resolve_mesh_path "${MESH_ETREE_RAW}" "${RUN_DIR_ABS}" "${BASE_IN_DIR}")"
    TMP_FIX_A="${TMP_PHASE_A_IN}.fix"
    TMP_FIX_B="${TMP_PHASE_B_IN}.fix"
    set_or_add_param "${TMP_PHASE_A_IN}" "${TMP_FIX_A}" "mesh_etree_output_file" "${EFFECTIVE_MESH_ETREE_PATH}"
    set_or_add_param "${TMP_PHASE_B_IN}" "${TMP_FIX_B}" "mesh_etree_output_file" "${EFFECTIVE_MESH_ETREE_PATH}"
    mv -f "${TMP_FIX_A}" "${TMP_PHASE_A_IN}"
    mv -f "${TMP_FIX_B}" "${TMP_PHASE_B_IN}"
  fi

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
echo "Run dir:        ${RUN_DIR}"
if [[ -n "${BASE_IN}" ]]; then
  echo "Base input:     ${BASE_IN}"
  echo "Phase A input:  ${PHASE_A_IN} (generated)"
  echo "Phase B input:  ${PHASE_B_IN} (generated)"
  if [[ -n "${EFFECTIVE_MESH_ETREE_PATH:-}" ]]; then
    echo "Mesh etree:     ${EFFECTIVE_MESH_ETREE_PATH}"
  fi
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

if [[ -n "${EFFECTIVE_MESH_ETREE_PATH:-}" ]]; then
  mkdir -p -- "$(dirname -- "${EFFECTIVE_MESH_ETREE_PATH}")"
fi

if [[ -n "${PHASE_A_NP}" ]]; then
  "${MPI_RUN}" -np "${PHASE_A_NP}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
else
  "${MPI_RUN}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
fi

if [[ -n "${EFFECTIVE_MESH_ETREE_PATH:-}" ]] && [[ ! -s "${EFFECTIVE_MESH_ETREE_PATH}" ]]; then
  echo "ERROR: Phase A did not produce a readable mesh etree file: ${EFFECTIVE_MESH_ETREE_PATH}" >&2
  exit 1
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
