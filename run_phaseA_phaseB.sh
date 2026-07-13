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
: "${HERCULES_SPARSE_CLEAR_MODE:=full}"
: "${HERCULES_FULL_CLEAR_OMP:=1}"
: "${HERCULES_STATION_IO_OMP:=1}"
: "${HERCULES_STATION_ROUTE_MODE:=interval}"
: "${HERCULES_PLANE_IO_OMP:=1}"
: "${HERCULES_CONTRIB_PACK_OMP:=1}"
: "${HERCULES_CONTRIB_UNPACK_OMP:=1}"
: "${HERCULES_CONTRIB_RECV_MODE:=waitsome}"
: "${HERCULES_SPARSE_UPLOAD_MODE:=dense}"
: "${HERCULES_PHASE_MODE:=two-phase}"
: "${HERCULES_UNIFIED_WRITE_ARCHIVE:=1}"
if [[ -z "${HERCULES_MESH_IO_MODE+x}" ]]; then
  if [[ "${HERCULES_PHASE_MODE}" == "unified" ]]; then
    HERCULES_MESH_IO_MODE=mpiio
  else
    HERCULES_MESH_IO_MODE=sharded
  fi
fi
export OMP_PROC_BIND OMP_PLACES
export HERCULES_ADJUST_ASSIGN_OMP HERCULES_ADJUST_DIST_OMP
export HERCULES_SPARSE_CLEAR_MODE HERCULES_FULL_CLEAR_OMP
export HERCULES_STATION_IO_OMP HERCULES_STATION_ROUTE_MODE
export HERCULES_PLANE_IO_OMP
export HERCULES_CONTRIB_PACK_OMP HERCULES_CONTRIB_UNPACK_OMP
export HERCULES_CONTRIB_RECV_MODE
export HERCULES_SPARSE_UPLOAD_MODE
export HERCULES_MESH_IO_MODE
export HERCULES_PHASE_MODE HERCULES_UNIFIED_WRITE_ARCHIVE

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
UNIFIED_IN=""
PHASE_A_MPI_ARGS="${PHASE_A_MPI_ARGS:-}"
PHASE_B_MPI_ARGS="${PHASE_B_MPI_ARGS:-}"

if [[ $# -ge 3 ]]; then
  BASE_IN="$3"
fi

case "${HERCULES_PHASE_MODE}" in
  two-phase|unified) ;;
  *)
    echo "ERROR: HERCULES_PHASE_MODE must be two-phase or unified" >&2
    exit 1
    ;;
esac
case "${HERCULES_UNIFIED_WRITE_ARCHIVE}" in
  0|1) ;;
  *)
    echo "ERROR: HERCULES_UNIFIED_WRITE_ARCHIVE must be 0 or 1" >&2
    exit 1
    ;;
esac
if [[ "${HERCULES_PHASE_MODE}" == "unified" && -z "${BASE_IN}" ]]; then
  echo "ERROR: unified mode requires BASE_IN (the third argument or BASE_IN env)" >&2
  exit 1
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
TMP_UNIFIED_IN=""

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

mpi_reports_cuda_support() {
  local ompi_info_bin="${OMPI_INFO:-ompi_info}"

  "${ompi_info_bin}" --parsable --all 2>/dev/null |
    grep -F 'mpi_built_with_cuda_support:value:true' >/dev/null
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
  TMP_PHASE_A_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phaseA.XXXXXX")"
  TMP_PHASE_B_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phaseB.XXXXXX")"
  TMP_UNIFIED_IN="$(mktemp "${TMPDIR_ROOT}/hercules_unified.XXXXXX")"

  cleanup_tmp_inputs() {
    [[ -n "${TMP_PHASE_A_IN}" ]] && rm -f "${TMP_PHASE_A_IN}"
    [[ -n "${TMP_PHASE_B_IN}" ]] && rm -f "${TMP_PHASE_B_IN}"
    [[ -n "${TMP_UNIFIED_IN}" ]] && rm -f "${TMP_UNIFIED_IN}"
  }
  trap cleanup_tmp_inputs EXIT

  make_phase_input "${BASE_IN}" "${TMP_PHASE_A_IN}" 1 1 0
  make_phase_input "${BASE_IN}" "${TMP_PHASE_B_IN}" 0 0 1
  make_phase_input "${BASE_IN}" "${TMP_UNIFIED_IN}" "${HERCULES_UNIFIED_WRITE_ARCHIVE}" 0 0

  # Normalize mesh etree path so both phases use the same absolute location.
  BASE_IN_DIR="$(cd -- "$(dirname -- "${BASE_IN}")" && pwd)"
  RUN_DIR_ABS="$(cd -- "${RUN_DIR}" && pwd)"
  MESH_ETREE_RAW="$(extract_config_value "${BASE_IN}" "mesh_etree_output_file")"
  MESH_ETREE_RAW="$(trim_ws "${MESH_ETREE_RAW}")"
  if [[ -z "${MESH_ETREE_RAW}" ]] &&
     { [[ "${HERCULES_PHASE_MODE}" == "two-phase" ]] ||
       [[ "${HERCULES_UNIFIED_WRITE_ARCHIVE}" == "1" ]]; }; then
    echo "ERROR: mesh_etree_output_file is required for persistent mesh output" >&2
    exit 1
  fi
  if [[ -n "${MESH_ETREE_RAW}" ]]; then
    EFFECTIVE_MESH_ETREE_PATH="$(resolve_mesh_path "${MESH_ETREE_RAW}" "${RUN_DIR_ABS}" "${BASE_IN_DIR}")"
    TMP_FIX_A="${TMP_PHASE_A_IN}.fix"
    TMP_FIX_B="${TMP_PHASE_B_IN}.fix"
    TMP_FIX_U="${TMP_UNIFIED_IN}.fix"
    set_or_add_param "${TMP_PHASE_A_IN}" "${TMP_FIX_A}" "mesh_etree_output_file" "${EFFECTIVE_MESH_ETREE_PATH}"
    set_or_add_param "${TMP_PHASE_B_IN}" "${TMP_FIX_B}" "mesh_etree_output_file" "${EFFECTIVE_MESH_ETREE_PATH}"
    set_or_add_param "${TMP_UNIFIED_IN}" "${TMP_FIX_U}" "mesh_etree_output_file" "${EFFECTIVE_MESH_ETREE_PATH}"
    mv -f "${TMP_FIX_A}" "${TMP_PHASE_A_IN}"
    mv -f "${TMP_FIX_B}" "${TMP_PHASE_B_IN}"
    mv -f "${TMP_FIX_U}" "${TMP_UNIFIED_IN}"
  fi

  PHASE_A_IN="${TMP_PHASE_A_IN}"
  PHASE_B_IN="${TMP_PHASE_B_IN}"
  UNIFIED_IN="${TMP_UNIFIED_IN}"
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

detect_alloc_cpus_per_node() {
  local raw=""

  if [[ -n "${SLURM_CPUS_ON_NODE:-}" ]] && [[ "${SLURM_CPUS_ON_NODE}" =~ ^[0-9]+$ ]]; then
    echo "${SLURM_CPUS_ON_NODE}"
    return
  fi
  if [[ -n "${SLURM_JOB_CPUS_PER_NODE:-}" ]]; then
    raw="${SLURM_JOB_CPUS_PER_NODE%%,*}"
    raw="${raw%%(*}"
    if [[ "${raw}" =~ ^[0-9]+$ ]]; then
      echo "${raw}"
      return
    fi
  fi
  if [[ -n "${PBS_NCPUS:-}" ]] && [[ "${PBS_NCPUS}" =~ ^[0-9]+$ ]]; then
    echo "${PBS_NCPUS}"
    return
  fi
  if command -v getconf >/dev/null 2>&1; then
    raw="$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)"
    if [[ "${raw}" =~ ^[0-9]+$ ]]; then
      echo "${raw}"
      return
    fi
  fi
  echo "0"
}

detect_local_numa_nodes() {
  local raw=""

  if command -v lscpu >/dev/null 2>&1; then
    raw="$(lscpu -p=NODE 2>/dev/null | awk -F, '!/^#/ && $1 ~ /^[0-9]+$/ { seen[$1] = 1 } END { print length(seen) }')"
    if [[ "${raw}" =~ ^[0-9]+$ ]] && [[ "${raw}" -gt 0 ]]; then
      echo "${raw}"
      return
    fi
  fi

  if command -v numactl >/dev/null 2>&1; then
    raw="$(numactl --hardware 2>/dev/null | awk '/available:[[:space:]]+[0-9]+[[:space:]]+nodes/ { print $2; exit }')"
    if [[ "${raw}" =~ ^[0-9]+$ ]] && [[ "${raw}" -gt 0 ]]; then
      echo "${raw}"
      return
    fi
  fi

  echo "1"
}

ALLOC_NODES="$(detect_alloc_nodes)"
ALLOC_CPUS_PER_NODE="$(detect_alloc_cpus_per_node)"
LOCAL_NUMA_NODES="$(detect_local_numa_nodes)"
ALLOC_TOTAL_CPUS=0

if [[ "${ALLOC_NODES}" -gt 0 ]] && [[ "${ALLOC_CPUS_PER_NODE}" -gt 0 ]]; then
  ALLOC_TOTAL_CPUS="$((ALLOC_NODES * ALLOC_CPUS_PER_NODE))"
fi

PHASE_A_OMP_THREADS="${PHASE_A_OMP_NUM_THREADS:-1}"
PHASE_B_OMP_THREADS="${PHASE_B_OMP_NUM_THREADS:-${OMP_NUM_THREADS:-}}"

if ! [[ "${PHASE_A_OMP_THREADS}" =~ ^[0-9]+$ ]] ||
   [[ "${PHASE_A_OMP_THREADS}" -lt 1 ]]; then
  echo "ERROR: PHASE_A_OMP_NUM_THREADS must be a positive integer" >&2
  exit 1
fi

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

if [[ "${ALLOC_NODES}" -gt 1 ]] && [[ "${PHASE_B_NP}" -gt "${ALLOC_NODES}" ]] &&
   (( PHASE_B_NP % ALLOC_NODES != 0 )); then
  echo "ERROR: Phase B requests ${PHASE_B_NP} MPI ranks on ${ALLOC_NODES} allocated nodes." >&2
  echo "       The rank count is not divisible by the node count, so the launcher cannot" >&2
  echo "       assign an even rank-per-node GPU layout." >&2
  echo "       Use PHASE_B_NP = allocated_nodes x GPUs_per_node, e.g. 40 ranks for" >&2
  echo "       10 GB200 nodes or 48 ranks for 12 GB200 nodes." >&2
  exit 1
fi

MPI_HAS_CUDA_SUPPORT=0
if mpi_reports_cuda_support; then
  MPI_HAS_CUDA_SUPPORT=1
fi

PHASE_A_UCX_TLS_VALUE=""
PHASE_B_UCX_TLS_VALUE=""
PHASE_A_UCX_TLS_SOURCE=""
PHASE_B_UCX_TLS_SOURCE=""
PHASE_B_OMPI_MCA_PML_VALUE=""
PHASE_B_OMPI_MCA_PML_SOURCE=""
PHASE_B_UNSET_OMPI_MCA_BTL=0
PHASE_B_OMPI_MCA_BTL_SOURCE=""

if [[ -n "${UCX_TLS+x}" ]]; then
  PHASE_A_UCX_TLS_VALUE="${UCX_TLS}"
  PHASE_B_UCX_TLS_VALUE="${UCX_TLS}"
  PHASE_A_UCX_TLS_SOURCE="user"
  PHASE_B_UCX_TLS_SOURCE="user"
else
  if [[ -n "${PHASE_A_UCX_TLS+x}" ]]; then
    PHASE_A_UCX_TLS_VALUE="${PHASE_A_UCX_TLS}"
    PHASE_A_UCX_TLS_SOURCE="phase-a user"
  elif [[ "${MPI_HAS_CUDA_SUPPORT}" -eq 1 ]]; then
    PHASE_A_UCX_TLS_VALUE="sm,self,tcp"
    PHASE_A_UCX_TLS_SOURCE="phase-a cpu default"
  fi

  if [[ -n "${PHASE_B_UCX_TLS+x}" ]]; then
    PHASE_B_UCX_TLS_VALUE="${PHASE_B_UCX_TLS}"
    PHASE_B_UCX_TLS_SOURCE="phase-b user"
  elif [[ "${MPI_HAS_CUDA_SUPPORT}" -eq 1 ]]; then
    PHASE_B_UCX_TLS_VALUE="dc,sm,self,cuda"
    PHASE_B_UCX_TLS_SOURCE="phase-b cuda-aware default"
  fi
fi

if [[ -n "${OMPI_MCA_pml+x}" ]]; then
  PHASE_B_OMPI_MCA_PML_VALUE="${OMPI_MCA_pml}"
  PHASE_B_OMPI_MCA_PML_SOURCE="user"
elif [[ -n "${PHASE_B_OMPI_MCA_PML+x}" ]]; then
  PHASE_B_OMPI_MCA_PML_VALUE="${PHASE_B_OMPI_MCA_PML}"
  PHASE_B_OMPI_MCA_PML_SOURCE="phase-b user"
elif [[ "${MPI_HAS_CUDA_SUPPORT}" -eq 1 ]]; then
  PHASE_B_OMPI_MCA_PML_VALUE="ucx"
  PHASE_B_OMPI_MCA_PML_SOURCE="phase-b cuda-aware default"
fi

if [[ "${MPI_HAS_CUDA_SUPPORT}" -eq 1 && "${PHASE_B_KEEP_OMPI_MCA_BTL:-0}" != "1" ]]; then
  PHASE_B_UNSET_OMPI_MCA_BTL=1
  PHASE_B_OMPI_MCA_BTL_SOURCE="phase-b cuda-aware default"
fi

if [[ "${ALLOC_TOTAL_CPUS}" -gt 0 ]] && [[ -n "${PHASE_A_NP}" ]] &&
   [[ "${PHASE_A_NP}" -gt "${ALLOC_TOTAL_CPUS}" ]]; then
  case " ${PHASE_A_MPI_ARGS:-} " in
    *" --map-by :OVERSUBSCRIBE "*|*" --oversubscribe "*)
      ;;
    *)
      echo "ERROR: Phase A requests ${PHASE_A_NP} MPI ranks, but the allocation provides only" >&2
      echo "       ${ALLOC_TOTAL_CPUS} CPU slots (${ALLOC_NODES} node(s) x ${ALLOC_CPUS_PER_NODE} CPUs/node)." >&2
      echo "       Phase A is CPU meshing, so on this allocation use PHASE_A_NP <= ${ALLOC_TOTAL_CPUS}," >&2
      echo "       request more nodes, or explicitly allow oversubscription with" >&2
      echo "       PHASE_A_MPI_ARGS=\"--map-by :OVERSUBSCRIBE\"." >&2
      exit 1
      ;;
  esac
fi

echo "=== Hercules phase launcher ==="
echo "MPI launcher:   ${MPI_RUN}"
if [[ -n "${PHASE_A_UCX_TLS_SOURCE}" ]]; then
  if [[ "${PHASE_A_UCX_TLS_SOURCE}" == "user" ]]; then
    echo "UCX TLS:        user (${UCX_TLS})"
  else
    echo "Phase A UCX TLS:${PHASE_A_UCX_TLS_SOURCE} (${PHASE_A_UCX_TLS_VALUE})"
    if [[ -n "${PHASE_B_UCX_TLS_SOURCE}" ]]; then
      echo "Phase B UCX TLS:${PHASE_B_UCX_TLS_SOURCE} (${PHASE_B_UCX_TLS_VALUE})"
    fi
  fi
fi
if [[ -n "${PHASE_B_OMPI_MCA_PML_SOURCE}" ]]; then
  echo "Phase B OMPI PML:${PHASE_B_OMPI_MCA_PML_SOURCE} (${PHASE_B_OMPI_MCA_PML_VALUE})"
fi
if [[ "${PHASE_B_UNSET_OMPI_MCA_BTL}" -eq 1 ]]; then
  echo "Phase B OMPI BTL:${PHASE_B_OMPI_MCA_BTL_SOURCE} (unset)"
fi
echo "psolve binary:  ${PSOLVE_BIN}"
echo "Run dir:        ${RUN_DIR}"
echo "Phase mode:     ${HERCULES_PHASE_MODE}"
echo "Mesh I/O mode:  ${HERCULES_MESH_IO_MODE}"
if [[ -n "${BASE_IN}" ]]; then
  echo "Base input:     ${BASE_IN}"
  echo "Phase A input:  ${PHASE_A_IN} (generated)"
  echo "Phase B input:  ${PHASE_B_IN} (generated)"
  if [[ "${HERCULES_PHASE_MODE}" == "unified" ]]; then
    echo "Unified input:  ${UNIFIED_IN} (generated)"
  fi
  if [[ -n "${EFFECTIVE_MESH_ETREE_PATH:-}" ]]; then
    echo "Mesh path:      ${EFFECTIVE_MESH_ETREE_PATH}"
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
if [[ "${ALLOC_CPUS_PER_NODE}" -gt 0 ]]; then
  echo "CPUs/node:      ${ALLOC_CPUS_PER_NODE}"
fi
if [[ -n "${PHASE_A_NP}" ]]; then
  echo "Phase A ranks:  ${PHASE_A_NP}"
else
  echo "Phase A ranks:  launcher default"
fi
echo "Phase B ranks:  ${PHASE_B_NP}"
if [[ -n "${PHASE_A_OMP_THREADS}" ]]; then
  echo "Phase A OMP:    ${PHASE_A_OMP_THREADS}"
fi
echo

if [[ "${HERCULES_PHASE_MODE}" == "two-phase" && -n "${PHASE_A_NP}" ]]; then
  echo "[Phase A] ${MPI_RUN} -np ${PHASE_A_NP} ${PSOLVE_BIN} ${PHASE_A_IN}"
elif [[ "${HERCULES_PHASE_MODE}" == "two-phase" ]]; then
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

if [[ -n "${PHASE_A_OMP_THREADS}" ]]; then
  A_ENV_VARS=("OMP_NUM_THREADS=${PHASE_A_OMP_THREADS}")
else
  A_ENV_VARS=()
fi
if [[ -n "${PHASE_A_UCX_TLS_SOURCE}" && "${PHASE_A_UCX_TLS_SOURCE}" != "user" ]]; then
  A_ENV_VARS+=("UCX_TLS=${PHASE_A_UCX_TLS_VALUE}")
fi
if [[ ${#A_ENV_VARS[@]} -gt 0 ]]; then
  A_ENV=(env "${A_ENV_VARS[@]}")
else
  A_ENV=()
fi

if [[ "${HERCULES_PHASE_MODE}" == "unified" ]]; then
  if [[ -z "${PHASE_A_NP}" ]]; then
    echo "ERROR: unified mode requires an explicit Phase A rank count" >&2
    exit 1
  fi
  if [[ -z "${PHASE_B_OMP_THREADS}" && "${ALLOC_NODES}" -gt 0 &&
        "${ALLOC_CPUS_PER_NODE}" -gt 0 ]] &&
     (( PHASE_B_NP % ALLOC_NODES == 0 )); then
    UNIFIED_RANKS_PER_NODE="$((PHASE_B_NP / ALLOC_NODES))"
    PHASE_B_OMP_THREADS="$((ALLOC_CPUS_PER_NODE / UNIFIED_RANKS_PER_NODE))"
  fi
  if [[ -z "${PHASE_B_OMP_THREADS}" ]] ||
     ! [[ "${PHASE_B_OMP_THREADS}" =~ ^[0-9]+$ ]] ||
     [[ "${PHASE_B_OMP_THREADS}" -lt 1 ]]; then
    echo "ERROR: set PHASE_B_OMP_NUM_THREADS to a positive integer for unified mode" >&2
    exit 1
  fi
  UNIFIED_MPI_ARGS_VALUE="${UNIFIED_MPI_ARGS:---bind-to none}"
  # shellcheck disable=SC2206
  U_ARGS=(${UNIFIED_MPI_ARGS_VALUE})
  U_ENV_VARS=(
    "OMP_NUM_THREADS=${PHASE_A_OMP_THREADS}"
    "OMP_PROC_BIND=false"
    "HERCULES_UNIFIED_PHASE_A_OMP_THREADS=${PHASE_A_OMP_THREADS}"
    "HERCULES_UNIFIED_PHASE_B_OMP_THREADS=${PHASE_B_OMP_THREADS}"
    "HERCULES_UNIFIED_PHASE_B_RANKS=${PHASE_B_NP}"
    "HERCULES_MESH_IO_MODE=${HERCULES_MESH_IO_MODE}"
  )
  if [[ "${PHASE_B_UNSET_OMPI_MCA_BTL}" -eq 1 ]]; then
    U_ENV_VARS=("-u" "OMPI_MCA_btl" "${U_ENV_VARS[@]}")
  fi
  if [[ -n "${PHASE_B_UCX_TLS_SOURCE}" &&
        "${PHASE_B_UCX_TLS_SOURCE}" != "user" ]]; then
    U_ENV_VARS+=("UCX_TLS=${PHASE_B_UCX_TLS_VALUE}")
  fi
  if [[ -n "${PHASE_B_OMPI_MCA_PML_SOURCE}" &&
        "${PHASE_B_OMPI_MCA_PML_SOURCE}" != "user" ]]; then
    U_ENV_VARS+=("OMPI_MCA_pml=${PHASE_B_OMPI_MCA_PML_VALUE}")
  fi
  echo "[Unified] ${MPI_RUN} -np ${PHASE_A_NP} ${U_ARGS[*]} ${PSOLVE_BIN} ${UNIFIED_IN}"
  env "${U_ENV_VARS[@]}" "${MPI_RUN}" -np "${PHASE_A_NP}" "${U_ARGS[@]}" \
    "${PSOLVE_BIN}" "${UNIFIED_IN}"
  if [[ "${HERCULES_UNIFIED_WRITE_ARCHIVE}" == "1" &&
        -n "${EFFECTIVE_MESH_ETREE_PATH:-}" &&
        ! -s "${EFFECTIVE_MESH_ETREE_PATH}" ]]; then
    echo "ERROR: unified run did not produce the MPI-IO mesh archive: ${EFFECTIVE_MESH_ETREE_PATH}" >&2
    exit 1
  fi
  echo
  echo "Unified mesh-and-solver run completed."
  exit 0
fi

if [[ -n "${PHASE_A_NP}" ]]; then
  "${A_ENV[@]}" "${MPI_RUN}" -np "${PHASE_A_NP}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
else
  "${A_ENV[@]}" "${MPI_RUN}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
fi

if [[ -n "${EFFECTIVE_MESH_ETREE_PATH:-}" ]] && [[ ! -s "${EFFECTIVE_MESH_ETREE_PATH}" ]]; then
  echo "ERROR: Phase A did not produce a readable mesh file or manifest: ${EFFECTIVE_MESH_ETREE_PATH}" >&2
  exit 1
fi

echo

PHASE_B_MAP_MODE="${PHASE_B_MAP_MODE:-auto}"
PHASE_B_MAP_ARGS=()
PHASE_B_MAP_DESC=""
PHASE_B_RANKS_PER_NODE=""
PHASE_B_RANKS_PER_NUMA=""

if [[ "${ALLOC_NODES}" -gt 0 ]] && [[ "${PHASE_B_NP}" -gt "${ALLOC_NODES}" ]]; then
  if (( PHASE_B_NP % ALLOC_NODES == 0 )); then
    PHASE_B_RANKS_PER_NODE="$((PHASE_B_NP / ALLOC_NODES))"
  fi
fi

if [[ -z "${PHASE_B_OMP_THREADS}" ]]; then
  if [[ -n "${PHASE_B_RANKS_PER_NODE}" ]] && [[ "${ALLOC_CPUS_PER_NODE}" -gt 0 ]]; then
    PHASE_B_OMP_THREADS="$((ALLOC_CPUS_PER_NODE / PHASE_B_RANKS_PER_NODE))"
  fi
fi

if [[ -n "${PHASE_B_RANKS_PER_NODE}" ]] && [[ "${ALLOC_CPUS_PER_NODE}" -gt 0 ]] &&
   [[ "${PHASE_B_OMP_THREADS:-}" =~ ^[0-9]+$ ]] && [[ "${PHASE_B_OMP_THREADS}" -gt 0 ]]; then
  PHASE_B_TOTAL_CPUS="$((PHASE_B_RANKS_PER_NODE * PHASE_B_OMP_THREADS))"
  if (( PHASE_B_TOTAL_CPUS > ALLOC_CPUS_PER_NODE )); then
    echo "ERROR: Phase B requests ${PHASE_B_RANKS_PER_NODE} ranks/node x ${PHASE_B_OMP_THREADS} OMP threads" >&2
    echo "       = ${PHASE_B_TOTAL_CPUS} CPUs/node, but only ${ALLOC_CPUS_PER_NODE} CPUs/node are available." >&2
    echo "       Set PHASE_B_OMP_NUM_THREADS to at most $((ALLOC_CPUS_PER_NODE / PHASE_B_RANKS_PER_NODE))." >&2
    exit 1
  fi
fi

case "${PHASE_B_MAP_MODE}" in
  per-node)
    if [[ "${ALLOC_NODES}" -gt 0 ]] && [[ "${PHASE_B_NP}" -gt "${ALLOC_NODES}" ]]; then
      echo "ERROR: Phase B map mode 'per-node' allows at most one rank per allocated node." >&2
      echo "       Requested Phase B ranks: ${PHASE_B_NP}" >&2
      echo "       Allocated nodes:         ${ALLOC_NODES}" >&2
      echo "       Use PHASE_B_MAP_MODE=auto for multi-GPU nodes." >&2
      exit 1
    fi
    PHASE_B_MAP_ARGS=(--map-by ppr:1:node --rank-by node)
    PHASE_B_MAP_DESC="per-node"
    ;;
  single-node)
    if [[ "${PHASE_B_OMP_THREADS:-}" =~ ^[0-9]+$ ]] && [[ "${PHASE_B_OMP_THREADS}" -gt 0 ]]; then
      PHASE_B_MAP_ARGS=(--map-by "slot:PE=${PHASE_B_OMP_THREADS}" --bind-to core --rank-by slot)
      PHASE_B_MAP_DESC="single-node slot binding (PE=${PHASE_B_OMP_THREADS})"
    else
      PHASE_B_MAP_ARGS=(--map-by slot --bind-to core --rank-by slot)
      PHASE_B_MAP_DESC="single-node slot binding"
    fi
    ;;
  numa)
    if [[ -z "${PHASE_B_RANKS_PER_NODE}" ]]; then
      echo "ERROR: Phase B map mode 'numa' requires a multi-rank-per-node Phase B layout." >&2
      echo "       Requested Phase B ranks: ${PHASE_B_NP}" >&2
      echo "       Allocated nodes:         ${ALLOC_NODES}" >&2
      exit 1
    fi
    if [[ "${LOCAL_NUMA_NODES}" -le 1 ]] ||
       (( PHASE_B_RANKS_PER_NODE % LOCAL_NUMA_NODES != 0 )); then
      echo "ERROR: Phase B map mode 'numa' requires ranks/node to divide evenly across local NUMA domains." >&2
      echo "       Ranks/node:        ${PHASE_B_RANKS_PER_NODE}" >&2
      echo "       Local NUMA domains:${LOCAL_NUMA_NODES}" >&2
      exit 1
    fi
    PHASE_B_RANKS_PER_NUMA="$((PHASE_B_RANKS_PER_NODE / LOCAL_NUMA_NODES))"
    if [[ "${PHASE_B_OMP_THREADS:-}" =~ ^[0-9]+$ ]] && [[ "${PHASE_B_OMP_THREADS}" -gt 0 ]]; then
      PHASE_B_MAP_ARGS=(--map-by "ppr:${PHASE_B_RANKS_PER_NUMA}:numa:PE=${PHASE_B_OMP_THREADS}" --bind-to core --rank-by slot)
      PHASE_B_MAP_DESC="numa(${PHASE_B_RANKS_PER_NUMA} ranks/numa x ${LOCAL_NUMA_NODES} numa, PE=${PHASE_B_OMP_THREADS})"
    else
      PHASE_B_MAP_ARGS=(--map-by "ppr:${PHASE_B_RANKS_PER_NUMA}:numa" --bind-to core --rank-by slot)
      PHASE_B_MAP_DESC="numa(${PHASE_B_RANKS_PER_NUMA} ranks/numa x ${LOCAL_NUMA_NODES} numa)"
    fi
    ;;
  auto)
    if [[ "${ALLOC_NODES}" -gt 1 ]] && [[ "${PHASE_B_NP}" -le "${ALLOC_NODES}" ]]; then
      PHASE_B_MAP_ARGS=(--map-by ppr:1:node --rank-by node)
      PHASE_B_MAP_DESC="auto(per-node)"
    elif [[ -n "${PHASE_B_RANKS_PER_NODE}" ]]; then
      if [[ "${PHASE_B_OMP_THREADS:-}" =~ ^[0-9]+$ ]] && [[ "${PHASE_B_OMP_THREADS}" -gt 0 ]]; then
        PHASE_B_MAP_ARGS=(--map-by "ppr:${PHASE_B_RANKS_PER_NODE}:node:PE=${PHASE_B_OMP_THREADS}" --bind-to core --rank-by slot)
        PHASE_B_MAP_DESC="auto(${PHASE_B_RANKS_PER_NODE} ranks/node, PE=${PHASE_B_OMP_THREADS})"
      else
        PHASE_B_MAP_ARGS=(--map-by "ppr:${PHASE_B_RANKS_PER_NODE}:node" --bind-to core --rank-by slot)
        PHASE_B_MAP_DESC="auto(${PHASE_B_RANKS_PER_NODE} ranks/node)"
      fi
    else
      if [[ "${PHASE_B_OMP_THREADS:-}" =~ ^[0-9]+$ ]] && [[ "${PHASE_B_OMP_THREADS}" -gt 0 ]]; then
        PHASE_B_MAP_ARGS=(--map-by "slot:PE=${PHASE_B_OMP_THREADS}" --bind-to core --rank-by slot)
        PHASE_B_MAP_DESC="auto(single-node, PE=${PHASE_B_OMP_THREADS})"
      else
        PHASE_B_MAP_ARGS=(--map-by slot --bind-to core --rank-by slot)
        PHASE_B_MAP_DESC="auto(single-node)"
      fi
    fi
    ;;
  *)
    echo "ERROR: invalid PHASE_B_MAP_MODE='${PHASE_B_MAP_MODE}'. Use auto|per-node|single-node|numa." >&2
    exit 1
    ;;
esac

if [[ -n "${PHASE_B_OMP_THREADS}" ]]; then
  echo "Phase B OMP:    ${PHASE_B_OMP_THREADS}"
fi
echo "Phase B mapping: ${PHASE_B_MAP_DESC}"
echo "[Phase B] ${MPI_RUN} -np ${PHASE_B_NP} ${PHASE_B_MAP_ARGS[*]} ${PSOLVE_BIN} ${PHASE_B_IN}"
if [[ -n "${PHASE_B_MPI_ARGS}" ]]; then
  # shellcheck disable=SC2206
  B_ARGS=(${PHASE_B_MPI_ARGS})
else
  B_ARGS=()
fi

B_ENV_VARS=()
if [[ "${PHASE_B_UNSET_OMPI_MCA_BTL}" -eq 1 ]]; then
  B_ENV_VARS+=("-u" "OMPI_MCA_btl")
fi
if [[ -n "${PHASE_B_OMP_THREADS}" ]]; then
  B_ENV_VARS+=("OMP_NUM_THREADS=${PHASE_B_OMP_THREADS}")
fi
if [[ -n "${PHASE_B_UCX_TLS_SOURCE}" &&
      "${PHASE_B_UCX_TLS_SOURCE}" != "user" ]]; then
  B_ENV_VARS+=("UCX_TLS=${PHASE_B_UCX_TLS_VALUE}")
fi
if [[ -n "${PHASE_B_OMPI_MCA_PML_SOURCE}" &&
      "${PHASE_B_OMPI_MCA_PML_SOURCE}" != "user" ]]; then
  B_ENV_VARS+=("OMPI_MCA_pml=${PHASE_B_OMPI_MCA_PML_VALUE}")
fi
if [[ ${#B_ENV_VARS[@]} -gt 0 ]]; then
  B_ENV=(env "${B_ENV_VARS[@]}")
else
  B_ENV=()
fi

"${B_ENV[@]}" "${MPI_RUN}" -np "${PHASE_B_NP}" "${PHASE_B_MAP_ARGS[@]}" "${B_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_B_IN}"
echo
echo "Two-phase run completed."
