#!/usr/bin/env bash
set -euo pipefail

samples_tsv="${1:-config/samples.tsv}"
env_file="${2:-config/pipeline.env}"

# This entry point preserves the original quality-filtering behavior and final
# output names. Only the IGH MaskPrimers stage is split into serial chunks.
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

# Chain-specific pRESTO options copied from the original workflow.
IGH_PRIMER_MAXLEN="${IGH_PRIMER_MAXLEN:-20}"
TRB_PRIMER_MAXLEN="${TRB_PRIMER_MAXLEN:-27}"
IGH_FILTER_NPROC="${IGH_FILTER_NPROC:-5}"
TRB_FILTER_NPROC="${TRB_FILTER_NPROC:-}"
IGH_MASK_NPROC="${IGH_MASK_NPROC:-1}"
TRB_MASK_NPROC="${TRB_MASK_NPROC:-}"
IGH_MASK_CHUNK_SIZE="${IGH_MASK_CHUNK_SIZE:-10000000}"

case "$IGH_MASK_CHUNK_SIZE" in
    ''|*[!0-9]*)
        echo "IGH_MASK_CHUNK_SIZE must be a positive integer." >&2
        exit 1
        ;;
esac
if [ "$IGH_MASK_CHUNK_SIZE" -lt 1 ]; then
    echo "IGH_MASK_CHUNK_SIZE must be a positive integer." >&2
    exit 1
fi

mkdir -p "$RESULTS_DIR" "$LOG_DIR"

count_fastq_reads() {
    local fastq="$1"
    awk '
        END {
            if (NR % 4 != 0) {
                print "Invalid FASTQ line count in " FILENAME ": " NR > "/dev/stderr"
                exit 1
            }
            print NR / 4
        }
    ' "$fastq"
}

run_mask_primers() {
    local sequence_file="$1"
    local primer_fasta="$2"
    local primer_maxlen="$3"
    local mask_nproc="$4"
    local log_file="$5"
    local outname="$6"

    local mask_args=(
        align
        -s "$sequence_file"
        -p "$primer_fasta"
        --skiprc
        --failed
        --mode tag
        --barcode
        --bf BARCODE
        --pf PRIMER
        --maxerror "$PRIMER_MAX_ERROR"
        --maxlen "$primer_maxlen"
        --log "$log_file"
        --outname "$outname"
    )
    if [ -n "$mask_nproc" ]; then
        mask_args+=(--nproc "$mask_nproc")
    fi
    MaskPrimers.py "${mask_args[@]}"
}

run_chunked_mask_primers() {
    local sample_id="$1"
    local sample_dir="$2"
    local read_label="$3"
    local quality_pass="$4"
    local primer_fasta="$5"
    local primer_maxlen="$6"
    local mask_nproc="$7"

    local work_dir="${sample_dir}/.${read_label}_mask_chunks"
    local input_dir="${work_dir}/input"
    local output_dir="${work_dir}/output"
    local merged_dir="${work_dir}/merged"
    local merged_pass="${merged_dir}/${read_label}_primers-pass.fastq"
    local merged_fail="${merged_dir}/${read_label}_primers-fail.fastq"
    local merged_log="${merged_dir}/${sample_id}_${read_label}.primer.log"

    if [ -e "$work_dir" ]; then
        echo "Chunk work directory already exists: ${work_dir}" >&2
        return 1
    fi

    mkdir -p "$input_dir" "$output_dir" "$merged_dir"
    SplitSeq.py count \
        -s "$quality_pass" \
        -n "$IGH_MASK_CHUNK_SIZE" \
        --outdir "$input_dir"

    local chunk_files=()
    mapfile -t chunk_files < <(
        find "$input_dir" -maxdepth 1 -type f -name '*_part*.fastq' -print | LC_ALL=C sort
    )
    if [ "${#chunk_files[@]}" -eq 0 ]; then
        echo "SplitSeq produced no FASTQ chunks for ${quality_pass}." >&2
        return 1
    fi

    : > "$merged_pass"
    : > "$merged_fail"
    : > "$merged_log"

    local original_reads
    original_reads=$(count_fastq_reads "$quality_pass")
    local total_input=0
    local total_pass=0
    local total_fail=0
    local chunk chunk_name part_id part_dir prefix
    local pass_file fail_file log_file required_file
    local input_reads pass_reads fail_reads

    for chunk in "${chunk_files[@]}"; do
        chunk_name=$(basename "$chunk" .fastq)
        part_id="${chunk_name##*_}"
        part_dir="${output_dir}/${part_id}"
        prefix="${part_dir}/${part_id}"
        pass_file="${prefix}_primers-pass.fastq"
        fail_file="${prefix}_primers-fail.fastq"
        log_file="${part_dir}/${part_id}.primer.log"
        mkdir -p "$part_dir"

        echo "Running MaskPrimers for ${sample_id} ${read_label} ${part_id}"
        run_mask_primers \
            "$chunk" "$primer_fasta" "$primer_maxlen" "$mask_nproc" \
            "$log_file" "$prefix"

        for required_file in "$pass_file" "$fail_file" "$log_file"; do
            if [ ! -f "$required_file" ]; then
                echo "Expected chunk output missing: ${required_file}" >&2
                return 1
            fi
        done

        input_reads=$(count_fastq_reads "$chunk")
        pass_reads=$(count_fastq_reads "$pass_file")
        fail_reads=$(count_fastq_reads "$fail_file")
        if [ $((pass_reads + fail_reads)) -ne "$input_reads" ]; then
            echo "Chunk count mismatch for ${part_id}: input=${input_reads}, pass=${pass_reads}, fail=${fail_reads}" >&2
            return 1
        fi

        cat "$pass_file" >> "$merged_pass"
        cat "$fail_file" >> "$merged_fail"
        cat "$log_file" >> "$merged_log"
        total_input=$((total_input + input_reads))
        total_pass=$((total_pass + pass_reads))
        total_fail=$((total_fail + fail_reads))
    done

    local merged_pass_reads merged_fail_reads
    merged_pass_reads=$(count_fastq_reads "$merged_pass")
    merged_fail_reads=$(count_fastq_reads "$merged_fail")
    if [ "$total_input" -ne "$original_reads" ] || \
       [ "$merged_pass_reads" -ne "$total_pass" ] || \
       [ "$merged_fail_reads" -ne "$total_fail" ] || \
       [ $((merged_pass_reads + merged_fail_reads)) -ne "$original_reads" ]; then
        echo "Merged count validation failed for ${sample_id} ${read_label}." >&2
        return 1
    fi

    mv -f -- "$merged_pass" "${sample_dir}/${read_label}_primers-pass.fastq"
    mv -f -- "$merged_fail" "${sample_dir}/${read_label}_primers-fail.fastq"
    mv -f -- "$merged_log" "${LOG_DIR}/${sample_id}_${read_label}.primer.log"

    local expected_work_dir="${sample_dir}/.${read_label}_mask_chunks"
    if [ "$work_dir" != "$expected_work_dir" ] || \
       [ "$work_dir" = "$sample_dir" ] || [ "$work_dir" = "/" ]; then
        echo "Refusing to remove unsafe work directory: ${work_dir}" >&2
        return 1
    fi
    rm -rf -- "$work_dir"
    echo "Validated and merged ${sample_id} ${read_label}: input=${original_reads}, pass=${merged_pass_reads}, fail=${merged_fail_reads}"
}

run_sample() {
    local sample_id="$1"
    local chain="$2"
    local read1="$3"
    local read2="$4"
    local primer_fasta="$5"

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

    if [ "$chain_upper" = "IGH" ] && [ "$MAX_PARALLEL_TASKS" -ne 1 ]; then
        echo "Chunked IGH processing requires MAX_PARALLEL_TASKS=1." >&2
        return 1
    fi

    for read_path in "$read1" "$read2"; do
        local read_label
        read_label=$(basename "$read_path")
        read_label="${read_label%.gz}"
        read_label="${read_label%.fastq}"
        read_label="${read_label%.fq}"

        # Refuse before FilterSeq can overwrite the quality-pass input that an
        # interrupted chunk directory was built from. The chunk function checks
        # again before creating the directory to guard against a later race.
        if [ "$chain_upper" = "IGH" ] && \
           [ -e "${sample_dir}/.${read_label}_mask_chunks" ]; then
            echo "Chunk work directory already exists: ${sample_dir}/.${read_label}_mask_chunks" >&2
            return 1
        fi

        # Step 1 is intentionally identical to the original script.
        local filter_args=(
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

        local quality_pass="${sample_dir}/${read_label}_quality-pass.fastq"
        if [ "$chain_upper" = "IGH" ]; then
            run_chunked_mask_primers \
                "$sample_id" "$sample_dir" "$read_label" "$quality_pass" \
                "$primer_fasta" "$primer_maxlen" "$mask_nproc"
        else
            run_mask_primers \
                "$quality_pass" "$primer_fasta" "$primer_maxlen" "$mask_nproc" \
                "${LOG_DIR}/${sample_id}_${read_label}.primer.log" \
                "${sample_dir}/${read_label}"
        fi
    done

    echo "Finished ${sample_id}"
}

wait_for_batch() {
    local status=0
    local pid
    for pid in "$@"; do
        if ! wait "$pid"; then
            status=1
        fi
    done
    return "$status"
}

task_pids=()
while IFS=$'\t' read -r sample_id chain read1 read2 primer_fasta mixcr_preset mixcr_mode; do
    [ -z "${sample_id:-}" ] && continue
    run_sample "$sample_id" "$chain" "$read1" "$read2" "$primer_fasta" &
    task_pids+=("$!")
    if [ "${#task_pids[@]}" -ge "$MAX_PARALLEL_TASKS" ]; then
        wait_for_batch "${task_pids[@]}"
        task_pids=()
    fi
done < <(tail -n +2 "$samples_tsv")

if [ "${#task_pids[@]}" -gt 0 ]; then
    wait_for_batch "${task_pids[@]}"
fi
echo "All primer filtering jobs completed."
