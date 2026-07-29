#!/usr/bin/env bash
set -euo pipefail

samples_tsv="${1:-config/samples.tsv}"
env_file="${2:-config/pipeline.env}"

# Each row in samples.tsv defines one sample and the primer FASTA to use.
# The env file keeps local paths and tunable parameters out of the code.
if [ -f "$env_file" ]; then
    # shellcheck disable=SC1090
    source "$env_file"
fi

RESULTS_DIR="${RESULTS_DIR:-results}"
LOG_DIR="${LOG_DIR:-logs}"
PREPROCESS_DIR="${PREPROCESS_DIR:-preprocessing}"
QUALITY_THRESHOLD="${QUALITY_THRESHOLD:-20}"
PRIMER_MAX_ERROR="${PRIMER_MAX_ERROR:-0.3}"
MAX_PARALLEL_TASKS="${MAX_PARALLEL_TASKS:-2}"

# Chain-specific pRESTO options.
# IGH and TRB use different primer lengths in the original workflow because
# the primer sets have different expected aligned lengths after allowing the
# mini-primer variants.
IGH_PRIMER_MAXLEN="${IGH_PRIMER_MAXLEN:-20}"
TRB_PRIMER_MAXLEN="${TRB_PRIMER_MAXLEN:-27}"
IGH_FILTER_NPROC="${IGH_FILTER_NPROC:-5}"
TRB_FILTER_NPROC="${TRB_FILTER_NPROC:-}"
IGH_MASK_NPROC="${IGH_MASK_NPROC:-1}"
TRB_MASK_NPROC="${TRB_MASK_NPROC:-}"

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

run_sample() {
    local sample_id="$1"
    local chain="$2"
    local read1="$3"
    local read2="$4"
    local primer_fasta="$5"

    # All intermediate FASTQ files for one sample are kept together so that
    # later steps can infer filenames from sample_id.
    local sample_dir="${RESULTS_DIR}/${sample_id}/${PREPROCESS_DIR}"
    mkdir -p "$sample_dir"

    echo "Processing ${sample_id} (${chain})"

    local chain_upper
    chain_upper=$(printf '%s' "$chain" | tr '[:lower:]' '[:upper:]')

    local primer_maxlen
    local filter_nproc
    local mask_nproc
    case "$chain_upper" in
        IGH)
            primer_maxlen="$IGH_PRIMER_MAXLEN"
            filter_nproc="$IGH_FILTER_NPROC"
            mask_nproc="$IGH_MASK_NPROC"
            ;;
        TRB|TRA)
            primer_maxlen="$TRB_PRIMER_MAXLEN"
            filter_nproc="$TRB_FILTER_NPROC"
            mask_nproc="$TRB_MASK_NPROC"
            ;;
        *)
            echo "Unsupported chain '${chain}' for sample '${sample_id}'. Expected IGH, TRB, or TRA." >&2
            return 1
            ;;
    esac

    for read_path in "$read1" "$read2"; do
        local read_label
        read_label=$(basename "$read_path")
        read_label="${read_label%.gz}"
        read_label="${read_label%.fastq}"
        read_label="${read_label%.fq}"

        # Step 1: quality filtering. Reads below QUALITY_THRESHOLD are written
        # to pRESTO's failed output; the passing reads feed primer masking.
        filter_args=(
            quality
            -s "$read_path"
            --failed
            -q "$QUALITY_THRESHOLD"
            --outname "${sample_dir}/${read_label}"
        )
        if [ -n "$filter_nproc" ]; then
            filter_args+=(--nproc "$filter_nproc")
        fi
        FilterSeq.py "${filter_args[@]}"

        # Step 2: primer masking and primer/barcode tag extraction.
        # --mode tag records the matched primer name in the FASTQ description.
        # --barcode stores the primer-derived barcode field used later to split
        # reads by Fprimer/Rprimer orientation.
        mask_args=(
            align
            -s "${sample_dir}/${read_label}_quality-pass.fastq"
            -p "$primer_fasta"
            --skiprc
            --failed
            --mode tag
            --barcode
            --bf BARCODE
            --pf PRIMER
            --maxerror "$PRIMER_MAX_ERROR"
            --maxlen "$primer_maxlen"
            --log "${LOG_DIR}/${sample_id}_${read_label}.primer.log"
            --outname "${sample_dir}/${read_label}"
        )
        if [ -n "$mask_nproc" ]; then
            mask_args+=(--nproc "$mask_nproc")
        fi
        MaskPrimers.py "${mask_args[@]}"
    done

    echo "Finished ${sample_id}"
}

task_count=0
tail -n +2 "$samples_tsv" | while IFS=$'\t' read -r sample_id chain read1 read2 primer_fasta mixcr_preset mixcr_mode; do
    [ -z "${sample_id:-}" ] && continue
    run_sample "$sample_id" "$chain" "$read1" "$read2" "$primer_fasta" &
    task_count=$((task_count + 1))
    if [ "$task_count" -ge "$MAX_PARALLEL_TASKS" ]; then
        wait
        task_count=0
    fi
done

wait
echo "All primer filtering jobs completed."
