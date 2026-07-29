#!/usr/bin/env bash
set -euo pipefail

samples_tsv="${1:-config/samples.tsv}"
env_file="${2:-config/pipeline.env}"

if [ -f "$env_file" ]; then
    # shellcheck disable=SC1090
    source "$env_file"
fi

RESULTS_DIR="${RESULTS_DIR:-results}"
PREPROCESS_DIR="${PREPROCESS_DIR:-preprocessing}"
MIXCR_OUTDIR="${MIXCR_OUTDIR:-mixcr}"
ALIGNMENT_CHUNKSIZE="${ALIGNMENT_CHUNKSIZE:-200000}"

filter_one_sample() {
    local sample_id="$1"
    local chain="$2"
    local chain_upper
    chain_upper=$(printf '%s' "$chain" | tr '[:lower:]' '[:upper:]')

    local mixcr_dir="${RESULTS_DIR}/${sample_id}/${PREPROCESS_DIR}/${MIXCR_OUTDIR}"
    local input_tsv="${mixcr_dir}/${sample_id}_alignments.tsv"
    local output_tsv="${mixcr_dir}/${sample_id}-filtered_alignments_poly8A.tsv"

    if [ ! -f "$input_tsv" ]; then
        echo "Missing MiXCR alignments file, skipping: $input_tsv" >&2
        return 0
    fi

    echo "Filtering ${sample_id} alignments for chain ${chain_upper}"

    SAMPLE_ID="$sample_id" \
    CHAIN="$chain_upper" \
    INPUT_TSV="$input_tsv" \
    OUTPUT_TSV="$output_tsv" \
    ALIGNMENT_CHUNKSIZE="$ALIGNMENT_CHUNKSIZE" \
    python3 - <<'PY'
import os

import pandas as pd


chain = os.environ["CHAIN"]
input_tsv = os.environ["INPUT_TSV"]
output_tsv = os.environ["OUTPUT_TSV"]
chunksize = int(os.environ.get("ALIGNMENT_CHUNKSIZE", "200000"))

first_chunk = True
kept_rows = 0
total_rows = 0

# Read in chunks because exported MiXCR alignments can be large.
for chunk in pd.read_csv(input_tsv, sep="\t", chunksize=chunksize):
    total_rows += len(chunk)
    clone_id = pd.to_numeric(chunk["cloneId"], errors="coerce")
    # Keep reads assigned to the expected receptor chain and remove unassigned
    # alignments (cloneId == -1). This creates the input expected by the
    # spatial barcode and UMI matching step.
    filtered = chunk[
        chunk["chains"].astype(str).str.startswith(chain, na=False)
        & (clone_id != -1)
    ]
    kept_rows += len(filtered)
    filtered.to_csv(
        output_tsv,
        sep="\t",
        index=False,
        mode="w" if first_chunk else "a",
        header=first_chunk,
    )
    first_chunk = False

if first_chunk:
    pd.DataFrame().to_csv(output_tsv, sep="\t", index=False)

print(f"Filtered alignments: kept {kept_rows} / {total_rows} rows -> {output_tsv}")
PY
}

tail -n +2 "$samples_tsv" | while IFS=$'\t' read -r sample_id chain read1 read2 primer_fasta mixcr_preset mixcr_mode; do
    [ -z "${sample_id:-}" ] && continue
    filter_one_sample "$sample_id" "$chain"
done
