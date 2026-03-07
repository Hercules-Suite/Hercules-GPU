#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CASE_DIR="${CASE_DIR:-$PWD}"

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
PHASE_A_PREP_SOURCE="${PHASE_A_PREP_SOURCE:-0}"
PHASE_A_SOURCE_REUSE="${PHASE_A_SOURCE_REUSE:-0}"
PHASE_B_SOURCE_REUSE="${PHASE_B_SOURCE_REUSE:-1}"
PHASE_SOURCE_PREP="${PHASE_SOURCE_PREP:-1}"
PHASE_SOURCE_PREP_MPI_ARGS="${PHASE_SOURCE_PREP_MPI_ARGS:-}"
PHASE_VERIFY_SOURCE_FILES="${PHASE_VERIFY_SOURCE_FILES:-1}"
PHASE_SOURCE_PREP_IF_MISSING="${PHASE_SOURCE_PREP_IF_MISSING:-0}"
PHASE_PIN_HOSTS="${PHASE_PIN_HOSTS:-1}"
PHASE_HOSTFILE="${PHASE_HOSTFILE:-}"

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
TMP_PHASE_PREP_IN=""
TMP_PHASE_HOSTFILE=""

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

extract_config_value() {
  local src="$1"
  local key="$2"

  awk -F'=' -v key="${key}" '
    $0 ~ "^[[:space:]]*" key "[[:space:]]*=" {
      val = $2;
      sub(/[[:space:]]*#.*/, "", val);
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", val);
      found = val;
    }
    END {
      if (found != "") {
        print found;
      }
    }
  ' "${src}"
}

verify_source_files() {
  local input_file="$1"
  local rank_count="$2"
  local src_dir_raw=""
  local src_dir=""
  local missing_count=0
  local present_count=0
  local active_count=0
  local missing_list=""
  local r f sz

  src_dir_raw="$(extract_config_value "${input_file}" "source_directory_output")"
  if [[ -z "${src_dir_raw}" ]]; then
    echo "WARNING: could not parse source_directory_output from ${input_file}; skipping source file verification" >&2
    return 0
  fi

  if [[ "${src_dir_raw}" = /* ]]; then
    src_dir="${src_dir_raw}"
  elif [[ -d "${src_dir_raw}" ]]; then
    src_dir="${src_dir_raw}"
  else
    src_dir="${CASE_DIR}/${src_dir_raw}"
  fi

  if [[ ! -d "${src_dir}" ]]; then
    echo "ERROR: source_directory_output does not exist: ${src_dir}" >&2
    return 1
  fi

  for ((r = 0; r < rank_count; r++)); do
    f="${src_dir}/force_process.${r}"
    if [[ -f "${f}" ]]; then
      present_count=$((present_count + 1))
      sz="$(wc -c < "${f}" | tr -d ' ')"
      if [[ "${sz}" =~ ^[0-9]+$ ]] && [[ "${sz}" -gt 4 ]]; then
        active_count=$((active_count + 1))
      fi
    else
      missing_count=$((missing_count + 1))
      if [[ -z "${missing_list}" ]]; then
        missing_list="${r}"
      elif [[ $(echo "${missing_list}" | awk -F',' '{print NF}') -lt 8 ]]; then
        missing_list="${missing_list},${r}"
      fi
    fi
  done

  if [[ "${missing_count}" -gt 0 ]]; then
    echo "ERROR: Source Prep produced ${present_count}/${rank_count} force files in ${src_dir}" >&2
    echo "       Missing ranks: ${missing_list}" >&2
    return 1
  fi

  echo "Source file verification: ${present_count}/${rank_count} files present (${active_count} with force payload, ${rank_count} total ranks)"
  return 0
}

source_files_complete() {
  local input_file="$1"
  local rank_count="$2"
  local src_dir_raw=""
  local src_dir=""
  local r f

  src_dir_raw="$(extract_config_value "${input_file}" "source_directory_output")"
  if [[ -z "${src_dir_raw}" ]]; then
    return 1
  fi

  if [[ "${src_dir_raw}" = /* ]]; then
    src_dir="${src_dir_raw}"
  elif [[ -d "${src_dir_raw}" ]]; then
    src_dir="${src_dir_raw}"
  else
    src_dir="${CASE_DIR}/${src_dir_raw}"
  fi

  [[ -d "${src_dir}" ]] || return 1

  for ((r = 0; r < rank_count; r++)); do
    f="${src_dir}/force_process.${r}"
    [[ -f "${f}" ]] || return 1
  done

  return 0
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
    [[ -n "${TMP_PHASE_PREP_IN}" ]] && rm -f "${TMP_PHASE_PREP_IN}"
    [[ -n "${TMP_PHASE_HOSTFILE}" ]] && rm -f "${TMP_PHASE_HOSTFILE}"
  }
  trap cleanup_tmp_inputs EXIT

  make_phase_input "${BASE_IN}" "${TMP_PHASE_A_IN}" 1 1 0
  make_phase_input "${BASE_IN}" "${TMP_PHASE_B_IN}" 0 0 1
  TMP_PHASE_PREP_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phasePrep.XXXXXX.in")"
  make_phase_input "${BASE_IN}" "${TMP_PHASE_PREP_IN}" 0 1 1
  PHASE_A_IN="${TMP_PHASE_A_IN}"
  PHASE_B_IN="${TMP_PHASE_B_IN}"
  PHASE_PREP_IN="${TMP_PHASE_PREP_IN}"
else
  TMPDIR_ROOT="${TMPDIR:-/tmp}"
  TMP_PHASE_PREP_IN="$(mktemp "${TMPDIR_ROOT}/hercules_phasePrep.XXXXXX.in")"

  cleanup_tmp_inputs() {
    [[ -n "${TMP_PHASE_A_IN}" ]] && rm -f "${TMP_PHASE_A_IN}"
    [[ -n "${TMP_PHASE_B_IN}" ]] && rm -f "${TMP_PHASE_B_IN}"
    [[ -n "${TMP_PHASE_PREP_IN}" ]] && rm -f "${TMP_PHASE_PREP_IN}"
    [[ -n "${TMP_PHASE_HOSTFILE}" ]] && rm -f "${TMP_PHASE_HOSTFILE}"
  }
  trap cleanup_tmp_inputs EXIT

  make_phase_input "${PHASE_B_IN}" "${TMP_PHASE_PREP_IN}" 0 1 1
  PHASE_PREP_IN="${TMP_PHASE_PREP_IN}"
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

HOSTFILE_ARGS=()
if [[ -n "${PHASE_HOSTFILE}" ]]; then
  if [[ ! -f "${PHASE_HOSTFILE}" ]]; then
    echo "ERROR: PHASE_HOSTFILE not found: ${PHASE_HOSTFILE}" >&2
    exit 1
  fi
  HOSTFILE_ARGS=(--hostfile "${PHASE_HOSTFILE}")
elif [[ "${PHASE_PIN_HOSTS}" == "1" ]]; then
  TMPDIR_ROOT="${TMPDIR:-/tmp}"
  TMP_PHASE_HOSTFILE="$(mktemp "${TMPDIR_ROOT}/hercules_phase_hosts.XXXXXX")"

  if [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]]; then
    awk '!seen[$0]++ { print $0 }' "${PBS_NODEFILE}" >"${TMP_PHASE_HOSTFILE}"
  elif [[ -n "${SLURM_NODELIST:-}" ]] && command -v scontrol >/dev/null 2>&1; then
    scontrol show hostnames "${SLURM_NODELIST}" >"${TMP_PHASE_HOSTFILE}"
  fi

  if [[ -s "${TMP_PHASE_HOSTFILE}" ]]; then
    HOSTFILE_ARGS=(--hostfile "${TMP_PHASE_HOSTFILE}")
  else
    rm -f "${TMP_PHASE_HOSTFILE}"
    TMP_PHASE_HOSTFILE=""
  fi
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
echo "Phase A source prep: HERCULES_MESH_ONLY_PREP_SOURCE=${PHASE_A_PREP_SOURCE}"
echo "Phase B prep source substep: ${PHASE_SOURCE_PREP} (uses Phase B ranks/layout)"
echo "Phase B prep source if missing: ${PHASE_SOURCE_PREP_IF_MISSING}"
echo "Source file verification: ${PHASE_VERIFY_SOURCE_FILES}"
echo "Phase A source reuse: HERCULES_REUSE_SOURCE_FILES=${PHASE_A_SOURCE_REUSE}"
echo "Phase B source reuse: HERCULES_REUSE_SOURCE_FILES=${PHASE_B_SOURCE_REUSE}"
if [[ ${#HOSTFILE_ARGS[@]} -gt 0 ]]; then
  echo "Pinned hostfile: ${HOSTFILE_ARGS[1]}"
else
  echo "Pinned hostfile: disabled"
fi
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
  HERCULES_MESH_ONLY_PREP_SOURCE="${PHASE_A_PREP_SOURCE}" \
  HERCULES_REUSE_SOURCE_FILES="${PHASE_A_SOURCE_REUSE}" \
  "${MPI_RUN}" -np "${PHASE_A_NP}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
else
  HERCULES_MESH_ONLY_PREP_SOURCE="${PHASE_A_PREP_SOURCE}" \
  HERCULES_REUSE_SOURCE_FILES="${PHASE_A_SOURCE_REUSE}" \
  "${MPI_RUN}" "${A_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_A_IN}"
fi

echo
echo "[Phase B] ${MPI_RUN} -np ${PHASE_B_NP} --map-by ppr:1:node --rank-by node ${HOSTFILE_ARGS[*]} ${PSOLVE_BIN} ${PHASE_B_IN}"
if [[ -n "${PHASE_B_MPI_ARGS}" ]]; then
  # shellcheck disable=SC2206
  B_ARGS=(${PHASE_B_MPI_ARGS})
else
  B_ARGS=()
fi

if [[ "${PHASE_SOURCE_PREP}" == "1" ]]; then
  skip_source_prep=0
  if [[ "${PHASE_SOURCE_PREP_IF_MISSING}" == "1" ]] && source_files_complete "${PHASE_PREP_IN}" "${PHASE_B_NP}"; then
    skip_source_prep=1
    echo
    echo "[Phase B Prep: Source] skipped (required force_process files already exist for ${PHASE_B_NP} ranks)"
  fi

  if [[ "${skip_source_prep}" == "0" ]]; then
    echo
    echo "[Phase B Prep: Source] ${MPI_RUN} -np ${PHASE_B_NP} --map-by ppr:1:node --rank-by node ${HOSTFILE_ARGS[*]} ${PSOLVE_BIN} ${PHASE_PREP_IN}"
    if [[ -n "${PHASE_SOURCE_PREP_MPI_ARGS}" ]]; then
      # shellcheck disable=SC2206
      P_ARGS=(${PHASE_SOURCE_PREP_MPI_ARGS})
    else
      P_ARGS=()
    fi

    HERCULES_MESH_ONLY_PREP_SOURCE=1 \
    HERCULES_REUSE_SOURCE_FILES=0 \
    "${MPI_RUN}" -np "${PHASE_B_NP}" --map-by ppr:1:node --rank-by node "${HOSTFILE_ARGS[@]}" "${P_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_PREP_IN}"
  fi

  if [[ "${PHASE_VERIFY_SOURCE_FILES}" == "1" ]]; then
    verify_source_files "${PHASE_PREP_IN}" "${PHASE_B_NP}"
  fi
fi

HERCULES_REUSE_SOURCE_FILES="${PHASE_B_SOURCE_REUSE}" \
"${MPI_RUN}" -np "${PHASE_B_NP}" --map-by ppr:1:node --rank-by node "${HOSTFILE_ARGS[@]}" "${B_ARGS[@]}" "${PSOLVE_BIN}" "${PHASE_B_IN}"

echo
echo "Two-phase run completed."
