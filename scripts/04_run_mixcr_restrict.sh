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
MIXCR_SPECIES="${MIXCR_SPECIES:-hsa}"
MIXCR_OUTDIR="${MIXCR_OUTDIR:-mixcr}"
MIXCR_PRESET="${MIXCR_PRESET:-generic-amplicon}"
MIXCR_THREADS="${MIXCR_THREADS:-20}"
MIXCR_CPUSET="${MIXCR_CPUSET:-0-19}"
MIXCR_HEAP_MB="${MIXCR_HEAP_MB:-50000}"
MIXCR_INITIAL_HEAP_MB="${MIXCR_INITIAL_HEAP_MB:-8000}"

mixcr_cmd=(
    taskset -c "$MIXCR_CPUSET"
    mixcr
    "-Xmx${MIXCR_HEAP_MB}m"
    "-Xms${MIXCR_INITIAL_HEAP_MB}m"
)

run_mixcr_sample() {
    local sample_id="$1"
    local sample_dir="${RESULTS_DIR}/${sample_id}/${PREPROCESS_DIR}"
    local mixcr_dir="${sample_dir}/${MIXCR_OUTDIR}"
    local input_fastq="${sample_dir}/${sample_id}_concatenated.fastq"
    local prefix="${mixcr_dir}/${sample_id}"

    mkdir -p "$mixcr_dir"
    echo "Running MiXCR for ${sample_id}"

    # MiXCR annotation.
    # - generic-amplicon is used for the concatenated amplicon-style reads.
    # - VEnd-to-C floating boundaries match the retained manuscript workflow.
    # - saveOriginalReads keeps read descriptions available in exportAlignments.
    "${mixcr_cmd[@]}" align \
        -p "$MIXCR_PRESET" \
        --threads "$MIXCR_THREADS" \
        --species "$MIXCR_SPECIES" \
        --rna \
        --floating-left-alignment-boundary VEnd \
        --floating-right-alignment-boundary C \
        -OsaveOriginalReads=true \
        --trimming-quality-threshold 0 \
        -f \
        "$input_fastq" \
        "${prefix}.vdjca" \
        --report "${prefix}.report.txt"

    # Assemble clonotypes while keeping alignments, because downstream barcode
    # matching uses per-read target sequence, quality, chain, and cloneId fields.
    "${mixcr_cmd[@]}" assemble \
        --write-alignments \
        "${prefix}.vdjca" \
        -f \
        "${prefix}.clna" \
        --report "${prefix}.assemble.report.txt"

    # Export only columns required by the spatial barcode/UMI workflow and
    # downstream VDJC summaries.
    "${mixcr_cmd[@]}" exportAlignments \
        -cloneId \
        -descrsR1 \
        -targetSequences \
        -targetQualities \
        -vGene \
        -dGene \
        -jGene \
        -cGene \
        -chains \
        -nFeature CDR3 \
        -f \
        "${prefix}.clna" \
        "${prefix}_alignments.tsv"

    echo "Finished MiXCR for ${sample_id}"
}

tail -n +2 "$samples_tsv" | while IFS=$'\t' read -r sample_id chain read1 read2 primer_fasta mixcr_preset mixcr_mode; do
    [ -z "${sample_id:-}" ] && continue
    run_mixcr_sample "$sample_id"
done
