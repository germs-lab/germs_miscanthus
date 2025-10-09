#!/usr/bin/env bash
set -Eeuo pipefail

# Defaults
IN_DIR=""
OUT_DIR=""
THREADS=4
BASENAME="all"
CAT_MODE=0
R1PAT="R1"
R2PAT="R2"

usage() {
  cat <<EOF
Usage: $(basename "$0") --in IN_DIR --out OUT_DIR [--threads N] [--basename NAME] [--r1-pattern R1] [--r2-pattern R2] [--cat]

Modes:
  - Default (no --cat): Run FastQC on each FASTQ file individually.
  - --cat: Concatenate all R1 files together and all R2 files together, then run FastQC on the two concatenated files.

Options:
  --in, -i           Input directory with FASTQ/FASTQ.GZ files.
  --out, -o          Output directory for concatenated files and/or FastQC reports.
  --threads, -t      Threads for FastQC (default: ${THREADS}).
  --basename         Base name for concatenated outputs (default: ${BASENAME}).
  --r1-pattern       Pattern to identify R1 files (default: ${R1PAT}).
  --r2-pattern       Pattern to identify R2 files (default: ${R2PAT}).
  --cat              Enable concatenation mode.
  --help, -h         Show this help.

Notes:
- Supports .fastq and .fastq.gz.
- If you hit the shell arg limit with many files, see the xargs example near the end.
EOF
}

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in|-i) IN_DIR="${2:-}"; shift ;;
    --out|-o) OUT_DIR="${2:-}"; shift ;;
    --threads|-t) THREADS="${2:-}"; shift ;;
    --basename) BASENAME="${2:-}"; shift ;;
    --r1-pattern) R1PAT="${2:-}"; shift ;;
    --r2-pattern) R2PAT="${2:-}"; shift ;;
    --cat) CAT_MODE=1 ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
  esac
  shift
done

if [[ -z "${IN_DIR}" || -z "${OUT_DIR}" ]]; then
  echo "Error: --in and --out are required." >&2
  usage
  exit 1
fi

mkdir -p "${OUT_DIR}"

# Collect files (sorted for determinism)
mapfile -t R1_FILES < <(find "${IN_DIR}" -type f \( -name "*${R1PAT}*.fastq" -o -name "*${R1PAT}*.fastq.gz" \) | sort -V)
mapfile -t R2_FILES < <(find "${IN_DIR}" -type f \( -name "*${R2PAT}*.fastq" -o -name "*${R2PAT}*.fastq.gz" \) | sort -V)
mapfile -t ALL_FILES < <(find "${IN_DIR}" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) | sort -V)

if (( ${#ALL_FILES[@]} == 0 )); then
  echo "No FASTQ files found in ${IN_DIR}" >&2
  exit 2
fi

echo "Found ${#R1_FILES[@]} R1 files and ${#R2_FILES[@]} R2 files."

if (( CAT_MODE == 1 )); then
  # Concatenate and FastQC the merged files
  mkdir -p "${OUT_DIR}/qc_concat_R1" "${OUT_DIR}/qc_concat_R2"

  OUT_R1="${OUT_DIR}/qc_concat_R1/${BASENAME}_R1.fastq"
  OUT_R2="${OUT_DIR}/qc_concat_R2/${BASENAME}_R2.fastq"

  # Concatenate R1
  : > "${OUT_R1}"
  for f in "${R1_FILES[@]}"; do
    if [[ "${f}" == *.gz ]]; then
      zcat -- "${f}" >> "${OUT_R1}"
    else
      cat -- "${f}" >> "${OUT_R1}"
    fi
  done
  echo "Concatenated R1 -> ${OUT_R1}"

  # Concatenate R2
  : > "${OUT_R2}"
  for f in "${R2_FILES[@]}"; do
    if [[ "${f}" == *.gz ]]; then
      zcat -- "${f}" >> "${OUT_R2}"
    else
      cat -- "${f}" >> "${OUT_R2}"
    fi
  done
  echo "Concatenated R2 -> ${OUT_R2}"

  echo "Running FastQC (concatenated)..."
  fastqc -t "${THREADS}" -o "${OUT_DIR}/qc_concat_R1" "${OUT_R1}"
  fastqc -t "${THREADS}" -o "${OUT_DIR}/qc_concat_R2" "${OUT_R2}"

  echo "Done. Reports:"
  echo " - ${OUT_DIR}/qc_concat_R1"
  echo " - ${OUT_DIR}/qc_concat_R2"
else
# Ensure output dirs exist
  mkdir -p "${OUT_DIR}/qc_individual_R1" "${OUT_DIR}/qc_individual_R2"

  echo "Running FastQC on R1 files (${#R1_FILES[@]}) (threads=${THREADS})..."
    if (( ${#R1_FILES[@]} > 0 )); then
        fastqc -t "${THREADS}" -o "${OUT_DIR}/qc_individual_R1" "${R1_FILES[@]}"
  echo "Done. Reports: ${OUT_DIR}/qc_individual_R1"
else
  echo "No R1 files found." >&2
fi

echo "Running FastQC on R2 files (${#R2_FILES[@]}) (threads=${THREADS})..."
if (( ${#R2_FILES[@]} > 0 )); then
  fastqc -t "${THREADS}" -o "${OUT_DIR}/qc_individual_R2" "${R2_FILES[@]}"
  echo "Done. Reports: ${OUT_DIR}/qc_individual_R2"
else
  echo "No R2 files found." >&2
fi

  # If you have too many files for one command, use xargs instead:
  # find "${IN_DIR}" -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 \
  #   | xargs -0 -n 1 -P "${THREADS}" fastqc -o "${OUT_DIR}/qc_individual"

  echo "Done. Reports: ${OUT_DIR}/qc_individual"
fi