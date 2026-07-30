#!/usr/bin/env bash
#
# Decompress one dated TRB/IGH directory and generate its samples TSV.
#
# Usage:
#   bash 00_decompress_fq.sh /path/to/20260724/TRB
#   bash 00_decompress_fq.sh /path/to/20260724/IGH
#
set -euo pipefail

die() {
    echo "Error: $*" >&2
    exit 1
}

if [ "$#" -ne 1 ]; then
    die "Usage: bash $(basename "$0") <TRB-or-IGH-directory>"
fi

requested_workdir="$1"
[ -d "$requested_workdir" ] ||
    die "Directory does not exist: ${requested_workdir}"
workdir=$(cd "$requested_workdir" && pwd -P)

chain=$(basename "$workdir")
case "$chain" in
    TRB)
        primer_fasta="primers/TRB_primer_all.fasta"
        mixcr_mode="tcr"
        ;;
    IGH)
        primer_fasta="primers/IGH_primer_all.fasta"
        mixcr_mode="igh"
        ;;
    *)
        die "Unsupported chain directory '${chain}'; expected TRB or IGH"
        ;;
esac

dataset_id=$(basename "$(dirname "$workdir")")
[ -n "$dataset_id" ] || die "Unable to derive dataset ID from ${workdir}"
case "$dataset_id" in
    *$'\t'*|*$'\n'*|*$'\r'*)
        die "Dataset ID contains a tab or newline: ${dataset_id}"
        ;;
esac

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
default_project_root=$(cd "${script_dir}/.." && pwd -P)
project_root_input="${SPACIR_ROOT:-$default_project_root}"
[ -d "$project_root_input" ] ||
    die "SpaCir root does not exist: ${project_root_input}"
project_root=$(cd "$project_root_input" && pwd -P)

config_dir="${project_root}/config"
[ -d "$config_dir" ] ||
    die "Config directory does not exist: ${config_dir}"
[ -f "${project_root}/${primer_fasta}" ] ||
    die "Primer file does not exist: ${project_root}/${primer_fasta}"

samples_tsv="${config_dir}/${dataset_id}_${chain}_samples.tsv"
[ ! -e "$samples_tsv" ] ||
    die "Samples TSV already exists: ${samples_tsv}"

shopt -s nullglob
input_files=("${workdir}"/*.fq.gz)
shopt -u nullglob
[ "${#input_files[@]}" -gt 0 ] ||
    die "No .fq.gz files found in ${workdir}"

declare -A read1_by_sample=()
declare -A read2_by_sample=()
declare -A samples=()

for source_path in "${input_files[@]}"; do
    [ -f "$source_path" ] || continue
    filename=$(basename "$source_path")
    sample_id=""
    mate=""

    case "$filename" in
        *-"${chain}"_1.fq.gz)
            sample_id=${filename%_1.fq.gz}
            mate=1
            ;;
        *-"${chain}"_2.fq.gz)
            sample_id=${filename%_2.fq.gz}
            mate=2
            ;;
        *)
            die "Unsupported ${chain} FASTQ filename: ${filename}"
            ;;
    esac

    [ -n "$sample_id" ] ||
        die "Empty sample ID in FASTQ filename: ${filename}"
    case "$sample_id" in
        *$'\t'*|*$'\n'*|*$'\r'*)
            die "Sample ID contains a tab or newline: ${filename}"
            ;;
    esac

    samples["$sample_id"]=1
    if [ "$mate" -eq 1 ]; then
        [ -z "${read1_by_sample[$sample_id]+present}" ] ||
            die "Duplicate R1 for sample '${sample_id}'"
        read1_by_sample["$sample_id"]="$source_path"
    else
        [ -z "${read2_by_sample[$sample_id]+present}" ] ||
            die "Duplicate R2 for sample '${sample_id}'"
        read2_by_sample["$sample_id"]="$source_path"
    fi
done

[ "${#samples[@]}" -gt 0 ] ||
    die "No supported ${chain} FASTQ pairs found in ${workdir}"

mapfile -t ordered_samples < <(
    printf '%s\n' "${!samples[@]}" | LC_ALL=C sort -V
)

for sample_id in "${ordered_samples[@]}"; do
    [ -n "${read1_by_sample[$sample_id]+present}" ] ||
        die "Missing R1 for sample '${sample_id}'"
    [ -n "${read2_by_sample[$sample_id]+present}" ] ||
        die "Missing R2 for sample '${sample_id}'"
done

fastq_dir="${workdir}/fastq"
temp_dir=$(mktemp -d "${workdir}/.decompress_fq.XXXXXX")
mkdir "${temp_dir}/fastq"

published_paths=()
created_fastq_dir=0
completed=0

cleanup() {
    status=$?
    trap - EXIT

    if [ "$completed" -eq 0 ]; then
        for published_path in "${published_paths[@]}"; do
            case "$published_path" in
                "${fastq_dir}/"*.fastq|"${samples_tsv}")
                    rm -f -- "$published_path"
                    ;;
            esac
        done
        if [ "$created_fastq_dir" -eq 1 ] && [ -d "$fastq_dir" ]; then
            rmdir -- "$fastq_dir" 2>/dev/null || true
        fi
    fi

    if [ -n "${temp_dir:-}" ] && [ -d "$temp_dir" ]; then
        case "$temp_dir" in
            "${workdir}/.decompress_fq."*)
                rm -rf -- "$temp_dir"
                ;;
        esac
    fi

    exit "$status"
}
trap cleanup EXIT

processed=0
skipped=0
declare -a staged_paths=()
declare -a target_paths=()

for sample_id in "${ordered_samples[@]}"; do
    output_stem=${sample_id//-/_}
    for mate in 1 2; do
        if [ "$mate" -eq 1 ]; then
            source_path=${read1_by_sample[$sample_id]}
        else
            source_path=${read2_by_sample[$sample_id]}
        fi

        target_path="${fastq_dir}/${output_stem}_R${mate}.fastq"
        if [ -e "$target_path" ]; then
            [ -f "$target_path" ] ||
                die "FASTQ target exists but is not a file: ${target_path}"
            echo "Skip existing: ${target_path}"
            skipped=$((skipped + 1))
            continue
        fi

        staged_path="${temp_dir}/fastq/${output_stem}_R${mate}.fastq"
        echo "Decompress: ${source_path} -> ${target_path}"
        gzip -dc -- "$source_path" > "$staged_path"
        staged_paths+=("$staged_path")
        target_paths+=("$target_path")
        processed=$((processed + 1))
    done
done

staged_tsv="${temp_dir}/${dataset_id}_${chain}_samples.tsv"
printf 'sample_id\tchain\tread1\tread2\tprimer_fasta\tmixcr_preset\tmixcr_mode\n' \
    > "$staged_tsv"

for sample_id in "${ordered_samples[@]}"; do
    output_stem=${sample_id//-/_}
    read1_path="${fastq_dir}/${output_stem}_R1.fastq"
    read2_path="${fastq_dir}/${output_stem}_R2.fastq"
    printf '%s\t%s\t%s\t%s\t%s\tgeneric-amplicon\t%s\n' \
        "$sample_id" "$chain" "$read1_path" "$read2_path" \
        "$primer_fasta" "$mixcr_mode" >> "$staged_tsv"
done

if [ ! -d "$fastq_dir" ]; then
    mkdir "$fastq_dir"
    created_fastq_dir=1
fi

for index in "${!staged_paths[@]}"; do
    mv -T -- "${staged_paths[$index]}" "${target_paths[$index]}"
    published_paths+=("${target_paths[$index]}")
done

for sample_id in "${ordered_samples[@]}"; do
    output_stem=${sample_id//-/_}
    [ -f "${fastq_dir}/${output_stem}_R1.fastq" ] ||
        die "Missing decompressed R1 for sample '${sample_id}'"
    [ -f "${fastq_dir}/${output_stem}_R2.fastq" ] ||
        die "Missing decompressed R2 for sample '${sample_id}'"
done

mv -T -- "$staged_tsv" "$samples_tsv"
published_paths+=("$samples_tsv")

rm -rf -- "$temp_dir"
temp_dir=""
completed=1

echo "Completed ${chain}: decompressed=${processed}, skipped=${skipped}"
echo "Samples TSV: ${samples_tsv}"
