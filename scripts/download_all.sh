#!/usr/bin/env bash
# download_all.sh
# Usage: ./download_all.sh "<d>"
# Example: ./download_all.sh "1_cell_decrease_2_sigma/"

set -u

REMOTE="mameen@smatter-login.syr.edu"
REMOTE_BASE="/home/mameen"

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <d>"
  exit 1
fi

D_RAW="$1"
# Normalize to ensure exactly one trailing slash
if [[ "$D_RAW" == */ ]]; then
  D="$D_RAW"
else
  D="${D_RAW}/"
fi

# Make top-level local directory
mkdir -p "data/${D}"

# Copy top-level files
scp "${REMOTE}:${REMOTE_BASE}/${D}stresses.txt"       "data/${D}" || echo "warn: missing stresses.txt for ${D}"
scp "${REMOTE}:${REMOTE_BASE}/${D}histogram_data.csv" "data/${D}" || echo "warn: missing histogram_data.csv for ${D}"

# Copy per-run files for 000..099
for i in $(seq -w 0 99); do
  LOCAL_DIR="data/${D}${i}/"
  REMOTE_DIR="${REMOTE_BASE}/${D}${i}/"
  mkdir -p "$LOCAL_DIR"

  for f in errors.txt q_values.txt; do
    scp "${REMOTE}:${REMOTE_DIR}${f}" "$LOCAL_DIR" || echo "warn: missing ${f} for ${D}${i}"
  done
done
