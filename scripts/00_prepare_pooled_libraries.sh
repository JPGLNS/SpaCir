#!/usr/bin/env bash
set -euo pipefail

readonly work_dir="$PWD"
readonly trb_output="${work_dir}/TRB"
readonly igh_output="${work_dir}/IGH"
readonly mapping_output="${work_dir}/library_number_mapping.txt"

die() {
    echo "Error: $*" >&2
    exit 1
}

for output in "$trb_output" "$igh_output" "$mapping_output"; do
    if [ -e "$output" ]; then
        die "Output already exists: $(basename "$output")"
    fi
done

declare -A read1_by_library=()
declare -A read2_by_library=()
declare -A libraries=()

shopt -s nullglob
input_files=(*.fq.gz *.fastq.gz)
shopt -u nullglob

[ "${#input_files[@]}" -gt 0 ] ||
    die "No .fq.gz or .fastq.gz files found in ${work_dir}"

for input in "${input_files[@]}"; do
    [ -f "$input" ] || continue
    library_id=""
    mate=""

    case "$input" in
        *_raw_1.fq.gz)    library_id=${input%_raw_1.fq.gz}; mate=1 ;;
        *_raw_2.fq.gz)    library_id=${input%_raw_2.fq.gz}; mate=2 ;;
        *_raw_1.fastq.gz) library_id=${input%_raw_1.fastq.gz}; mate=1 ;;
        *_raw_2.fastq.gz) library_id=${input%_raw_2.fastq.gz}; mate=2 ;;
        *_1.fq.gz)        library_id=${input%_1.fq.gz}; mate=1 ;;
        *_2.fq.gz)        library_id=${input%_2.fq.gz}; mate=2 ;;
        *_1.fastq.gz)     library_id=${input%_1.fastq.gz}; mate=1 ;;
        *_2.fastq.gz)     library_id=${input%_2.fastq.gz}; mate=2 ;;
        *) die "Unsupported FASTQ filename: ${input}" ;;
    esac

    [ -n "$library_id" ] ||
        die "Empty library ID in FASTQ filename: ${input}"
    case "$library_id" in
        *$'\t'*|*$'\n'*|*$'\r'*)
            die "Library ID contains a tab or newline: ${input}"
            ;;
    esac

    libraries["$library_id"]=1
    if [ "$mate" -eq 1 ]; then
        [ -z "${read1_by_library[$library_id]+present}" ] ||
            die "Duplicate R1 files for library '${library_id}'"
        read1_by_library["$library_id"]="$input"
    else
        [ -z "${read2_by_library[$library_id]+present}" ] ||
            die "Duplicate R2 files for library '${library_id}'"
        read2_by_library["$library_id"]="$input"
    fi
done

[ "${#libraries[@]}" -gt 0 ] ||
    die "No supported paired FASTQ filenames found in ${work_dir}"

mapfile -t ordered_libraries < <(
    printf '%s\n' "${!libraries[@]}" | LC_ALL=C sort -V
)

for library_id in "${ordered_libraries[@]}"; do
    [ -n "${read1_by_library[$library_id]+present}" ] ||
        die "Missing R1 for library '${library_id}'"
    [ -n "${read2_by_library[$library_id]+present}" ] ||
        die "Missing R2 for library '${library_id}'"
done

temp_dir=$(mktemp -d "${work_dir}/.prepare_pooled_libraries.XXXXXX")
published_paths=()
completed=0

cleanup() {
    local status=$?
    trap - EXIT

    if [ "$completed" -eq 0 ]; then
        local path
        for path in "${published_paths[@]}"; do
            case "$path" in
                "$trb_output"|"$igh_output")
                    rm -rf -- "$path"
                    ;;
                "$mapping_output")
                    rm -f -- "$path"
                    ;;
            esac
        done
    fi

    if [ -n "${temp_dir:-}" ] && [ -d "$temp_dir" ]; then
        case "$temp_dir" in
            "${work_dir}"/.prepare_pooled_libraries.*)
                rm -rf -- "$temp_dir"
                ;;
        esac
    fi

    exit "$status"
}
trap cleanup EXIT

mkdir "${temp_dir}/TRB" "${temp_dir}/IGH"
printf 'library_id\tsample_number\n' > "${temp_dir}/library_number_mapping.txt"

sample_number=0
for library_id in "${ordered_libraries[@]}"; do
    sample_number=$((sample_number + 1))
    read1=${read1_by_library[$library_id]}
    read2=${read2_by_library[$library_id]}

    cp -- "$read1" "${temp_dir}/TRB/${sample_number}-TRB_1.fq.gz"
    cp -- "$read2" "${temp_dir}/TRB/${sample_number}-TRB_2.fq.gz"
    cp -- "$read1" "${temp_dir}/IGH/${sample_number}-IGH_1.fq.gz"
    cp -- "$read2" "${temp_dir}/IGH/${sample_number}-IGH_2.fq.gz"
    printf '%s\t%d\n' "$library_id" "$sample_number" \
        >> "${temp_dir}/library_number_mapping.txt"

    echo "Prepared ${library_id} -> ${sample_number}-TRB and ${sample_number}-IGH"
done

mv -T -- "${temp_dir}/TRB" "$trb_output"
published_paths+=("$trb_output")
mv -T -- "${temp_dir}/IGH" "$igh_output"
published_paths+=("$igh_output")
mv -T -- "${temp_dir}/library_number_mapping.txt" "$mapping_output"
published_paths+=("$mapping_output")

rmdir -- "$temp_dir"
temp_dir=""
completed=1

echo "Completed: ${#ordered_libraries[@]} libraries"
echo "Mapping: ${mapping_output}"
