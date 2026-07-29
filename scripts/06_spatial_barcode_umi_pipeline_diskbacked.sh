#!/usr/bin/env bash
set -euo pipefail

MAX_WORKERS="${MAX_WORKERS:-1}"
if [ "$MAX_WORKERS" != "1" ]; then
    echo "06_spatial_barcode_umi_pipeline_diskbacked.sh requires MAX_WORKERS=1" >&2
    exit 2
fi

STAGE6_CSV_FIELD_LIMIT="${STAGE6_CSV_FIELD_LIMIT:-268435456}"
if ! [[ "$STAGE6_CSV_FIELD_LIMIT" =~ ^[1-9][0-9]*$ ]]; then
    echo "STAGE6_CSV_FIELD_LIMIT must be a positive integer" >&2
    exit 2
fi

export MAX_WORKERS
export STAGE6_CSV_FIELD_LIMIT

BASE_DIR="${BASE_DIR:-results}"
SAMPLE_GLOB="${SAMPLE_GLOB:-*-TRB}"
PREPROCESS_DIR="${PREPROCESS_DIR:-preprocessing}"
MIXCR_OUTDIR="${MIXCR_OUTDIR:-mixcr}"

if [ ! -d "$BASE_DIR" ]; then
    echo "Stage 6 base directory does not exist: $BASE_DIR" >&2
    echo "Stage 6 sample summary: discovered=0 processed=0 failed=0"
    exit 1
fi

sample_manifest=$(mktemp) || {
    echo "Unable to create stage 6 sample discovery manifest" >&2
    exit 1
}
cleanup_sample_manifest() {
    rm -f -- "$sample_manifest"
}
trap cleanup_sample_manifest EXIT

if ! find "$BASE_DIR" -mindepth 1 -maxdepth 1 -type d \
    -name "$SAMPLE_GLOB" | sort > "$sample_manifest"; then
    echo "Stage 6 sample discovery failed: base=$BASE_DIR glob=$SAMPLE_GLOB" >&2
    echo "Stage 6 sample summary: discovered=0 processed=0 failed=0"
    exit 1
fi

sample_dirs=()
while IFS= read -r sample_dir; do
    sample_dirs+=("$sample_dir")
done < "$sample_manifest"
cleanup_sample_manifest
trap - EXIT

discovered_count=${#sample_dirs[@]}
processed_count=0
failed_count=0

if [ "$discovered_count" -eq 0 ]; then
    echo "No stage 6 samples matched: base=$BASE_DIR glob=$SAMPLE_GLOB" >&2
    echo "Stage 6 sample summary: discovered=0 processed=0 failed=0"
    exit 1
fi

for sample_dir in "${sample_dirs[@]}"; do
    sample_name=$(basename "$sample_dir")
    mixcr_dir="${sample_dir}/${PREPROCESS_DIR}/${MIXCR_OUTDIR}"
    echo "========================================"
    echo "Processing ${sample_name}"
    echo "Dir: ${mixcr_dir}"
    if [ ! -d "$mixcr_dir" ]; then
        echo "Missing MiXCR directory: sample=$sample_name path=$mixcr_dir" >&2
        failed_count=$((failed_count + 1))
        continue
    fi

    if (
    if ! cd "$mixcr_dir"; then
        echo "Unable to enter MiXCR directory: sample=$sample_name path=$mixcr_dir" >&2
        exit 1
    fi
SAMPLE_NAME="$sample_name" \
MIXCR_DIR="$PWD" \
python3 -u - <<'PY'
"""Disk-backed primitives for the bounded-memory stage 6 pipeline."""
import contextlib
import csv
import json
import math
import os
import re
import sqlite3
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd


DNA_BASES = "ACGT"
max_levenshtein_distances = [1]
umi_thresholds = [3]


class CounterState(dict):
    """Weak-referenceable counter map scoped to one output combination."""


def split_braced(value):
    if not isinstance(value, str) or value == "":
        return []
    return value.strip("{}").split("}{")


def join_braced(values):
    return "{" + "}{".join(map(str, values)) + "}"


def update_matrix_counters(row, cdr3_counts, umi_sets):
    """Accumulate final CDR3 observations and unique UMI support."""
    barcode = row.get("spatial_barcode_ori")
    cdr3_value = row.get("aaSeqCDR3_most")
    if not isinstance(barcode, str) or not isinstance(cdr3_value, str):
        return

    cdr3_values = split_braced(cdr3_value)
    barcode_counts = cdr3_counts.setdefault(barcode, {})
    for cdr3 in cdr3_values:
        barcode_counts[cdr3] = barcode_counts.get(cdr3, 0) + 1

    umi_value = row.get("UMI")
    if not isinstance(umi_value, str):
        return
    barcode_umis = umi_sets.setdefault(barcode, {})
    for cdr3, umi in zip(cdr3_values, split_braced(umi_value)):
        barcode_umis.setdefault(cdr3, set()).add(umi)


def calculate_shannon_from_counts(counts):
    try:
        from skbio.diversity.alpha import shannon
        return shannon(counts)
    except ImportError:
        values = [float(value) for value in counts if float(value) > 0]
        total = sum(values)
        return -sum(
            (value / total) * math.log(value / total)
            for value in values
        ) if total else 0.0


def update_gene_frequency_counters(row, gene_umi_sets, vj_umi_sets):
    """Accumulate indexed-compatible unique UMI support by gene and V-J pair."""
    umi_value = row.get("UMI")
    if not isinstance(umi_value, str):
        return
    umis = split_braced(umi_value)
    category_fields = (
        ("VGene", "bestVGene"),
        ("DGene", "bestDGene"),
        ("JGene", "bestJGene"),
        ("CGene", "bestCGene"),
    )
    for category, fieldname in category_fields:
        value = row.get(fieldname)
        if not isinstance(value, str):
            continue
        category_sets = gene_umi_sets.setdefault(category, {})
        for gene, umi in zip(split_braced(value), umis):
            category_sets.setdefault(gene, set()).add(umi)

    v_value = row.get("bestVGene")
    j_value = row.get("bestJGene")
    if not isinstance(v_value, str) or not isinstance(j_value, str):
        return
    v_genes = split_braced(v_value)
    j_genes = split_braced(j_value)
    if len(v_genes) == len(j_genes) == len(umis):
        for v_gene, j_gene, umi in zip(v_genes, j_genes, umis):
            vj_umi_sets.setdefault(f"{v_gene}_{j_gene}", set()).add(umi)


def calculate_mutation_rate_per_segment(segment):
    try:
        parts = segment.split("|")
        start, end = int(parts[0]), int(parts[1])
        sequence_length = end - start + 1
        if len(parts) > 5:
            mutation_count = len(re.findall(r"\d+", parts[5]))
            return (
                f"{mutation_count / sequence_length:.4f}",
                sequence_length,
                mutation_count,
            )
    except Exception:
        return "nan", 0, 0


def calculate_mutation_rate_and_details(value):
    if not isinstance(value, str) or value.strip("{}").lower() == "nan":
        return "{nan}", "{0}", "{0}"

    mutation_rates = []
    sequence_lengths = []
    mutation_counts = []
    for segment in split_braced(value):
        if segment.lower() == "nan":
            mutation_rates.append("nan")
            sequence_lengths.append("0")
            mutation_counts.append("0")
        else:
            rate, length, count = calculate_mutation_rate_per_segment(segment)
            mutation_rates.append(rate)
            sequence_lengths.append(str(length))
            mutation_counts.append(str(count))
    return (
        join_braced(mutation_rates),
        join_braced(sequence_lengths),
        join_braced(mutation_counts),
    )


def _calculate_region_rate(output_row, regions):
    mutation_count = sum(
        int(output_row[f"all{region}Alignments_mutation_count"])
        for region in regions
    )
    total_length = sum(
        int(output_row[f"all{region}Alignments_sequence_length"])
        for region in regions
    )
    return mutation_count / total_length if total_length > 0 else 0


def write_mutation_rows(
    row,
    summary_writer,
    expanded_writer,
    sample="<unknown>",
    stage="[8/8] mutation",
    input_file="<unknown>",
):
    """Write one barcode mutation summary and its expanded records immediately."""
    summary_row = dict(row)
    required_fields = (
        "spatial_barcode_ori",
        "aaSeqCDR3_most",
        "nSeqCDR3",
        "bestVGene",
        "bestDGene",
        "bestJGene",
        "bestCGene",
        "UMI",
        "allVAlignments",
        "allDAlignments",
        "allJAlignments",
        "allCAlignments",
    )
    for fieldname in required_fields:
        if fieldname not in summary_row:
            raise ValueError(
                f"{_row_error_context(row, sample, stage, fieldname, input_file)}; "
                "missing required mutation field"
            )

    for column in (
        "allVAlignments",
        "allDAlignments",
        "allJAlignments",
        "allCAlignments",
    ):
        rates, lengths, counts = calculate_mutation_rate_and_details(
            summary_row.get(column)
        )
        summary_row[f"{column}_mutation_rate"] = rates
        summary_row[f"{column}_sequence_length"] = lengths
        summary_row[f"{column}_mutation_count"] = counts

    _validate_parallel_braced_fields(
        summary_row,
        "aaSeqCDR3_most",
        sample,
        stage,
        input_file,
    )

    values_by_field = {
        "aaSeqCDR3_most": split_braced(summary_row["aaSeqCDR3_most"]),
        "nSeqCDR3": split_braced(summary_row["nSeqCDR3"]),
        "bestVGene": split_braced(summary_row["bestVGene"]),
        "bestDGene": split_braced(summary_row["bestDGene"]),
        "bestJGene": split_braced(summary_row["bestJGene"]),
        "bestCGene": split_braced(summary_row["bestCGene"]),
        "UMI": split_braced(summary_row["UMI"]),
        "allVAlignments_mutation_rate": split_braced(
            summary_row["allVAlignments_mutation_rate"]
        ),
        "allDAlignments_mutation_rate": split_braced(
            summary_row["allDAlignments_mutation_rate"]
        ),
        "allJAlignments_mutation_rate": split_braced(
            summary_row["allJAlignments_mutation_rate"]
        ),
        "allCAlignments_mutation_rate": split_braced(
            summary_row["allCAlignments_mutation_rate"]
        ),
        "allVAlignments_sequence_length": split_braced(
            summary_row["allVAlignments_sequence_length"]
        ),
        "allDAlignments_sequence_length": split_braced(
            summary_row["allDAlignments_sequence_length"]
        ),
        "allJAlignments_sequence_length": split_braced(
            summary_row["allJAlignments_sequence_length"]
        ),
        "allCAlignments_sequence_length": split_braced(
            summary_row["allCAlignments_sequence_length"]
        ),
        "allVAlignments_mutation_count": split_braced(
            summary_row["allVAlignments_mutation_count"]
        ),
        "allDAlignments_mutation_count": split_braced(
            summary_row["allDAlignments_mutation_count"]
        ),
        "allJAlignments_mutation_count": split_braced(
            summary_row["allJAlignments_mutation_count"]
        ),
        "allCAlignments_mutation_count": split_braced(
            summary_row["allCAlignments_mutation_count"]
        ),
    }
    summary_writer.writerow(summary_row)

    for index in range(len(values_by_field["aaSeqCDR3_most"])):
        output_row = {
            "spatial_barcode_ori": summary_row["spatial_barcode_ori"],
            **{
                fieldname: values[index]
                for fieldname, values in values_by_field.items()
            },
        }

        output_row["VDJ_mutation_rate"] = _calculate_region_rate(
            output_row, ("V", "D", "J")
        )
        output_row["VDJC_mutation_rate"] = _calculate_region_rate(
            output_row, ("V", "D", "J", "C")
        )
        output_row["VJ_mutation_rate"] = _calculate_region_rate(
            output_row, ("V", "J")
        )
        expanded_writer.writerow(output_row)


def extract_unique_barcode(spatial_barcode_ori):
    barcodes = re.findall(r"\{(.*?)\}", spatial_barcode_ori)
    return barcodes[0] if barcodes else None


def iter_csv_rows(path):
    """Yield one logical CSV row at a time, including quoted multiline cells."""
    field_limit = int(os.environ.get("STAGE6_CSV_FIELD_LIMIT", "268435456"))
    if field_limit <= 0:
        raise ValueError("STAGE6_CSV_FIELD_LIMIT must be a positive integer")
    csv.field_size_limit(field_limit)
    with Path(path).open(newline="", encoding="utf-8") as input_handle:
        yield from csv.DictReader(input_handle)


def _csv_fieldnames(path):
    with Path(path).open(newline="", encoding="utf-8") as input_handle:
        return list(csv.DictReader(input_handle).fieldnames or [])


def _pipeline_base(input_file):
    return str(input_file).split(".")[0]


def _sample_label(input_file):
    return Path(input_file).name.split(".", 1)[0]


def _append_fieldnames(fieldnames, additions):
    output = list(fieldnames)
    for fieldname in additions:
        if fieldname not in output:
            output.append(fieldname)
    return output


def _write_csv_rows_atomic(
    output_file,
    fieldnames,
    rows,
    write_header_if_empty=True,
):
    with atomic_output(output_file) as temporary_path:
        with temporary_path.open("w", newline="", encoding="utf-8") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=fieldnames,
                lineterminator=os.linesep,
            )
            row_iterator = iter(rows)
            wrote_header = False
            try:
                while True:
                    try:
                        row = next(row_iterator)
                    except StopIteration:
                        break
                    if row is not None:
                        if not wrote_header:
                            writer.writeheader()
                            wrote_header = True
                        writer.writerow(row)
                    del row
            finally:
                close_iterator = getattr(row_iterator, "close", None)
                if close_iterator is not None:
                    close_iterator()

            if not wrote_header:
                if write_header_if_empty:
                    writer.writeheader()
                else:
                    output_handle.write(os.linesep)


def _write_frame_atomic(frame, output_file, index=False, encoding=None):
    with atomic_output(output_file) as temporary_path:
        frame.to_csv(
            temporary_path,
            index=index,
            encoding=encoding,
            lineterminator=os.linesep,
        )


def _row_error_context(row, sample, stage, fieldname, input_file):
    barcode_value = row.get("spatial_barcode_ori", "")
    barcode = (
        row.get("spatial_bc")
        or extract_unique_barcode(barcode_value)
        or barcode_value
    )
    return (
        f"sample={sample}; stage={stage}; barcode={barcode}; "
        f"field={fieldname}; input_file={input_file}"
    )


def _validate_parallel_braced_fields(
    row,
    reference_field,
    sample,
    stage,
    input_file,
):
    reference_values = split_braced(row[reference_field])
    expected_length = len(reference_values)

    for fieldname, value in row.items():
        if fieldname == reference_field:
            continue
        if isinstance(value, str) and "{" in value and "}" in value:
            actual_length = len(split_braced(value))
            if actual_length != expected_length:
                raise ValueError(
                    f"{_row_error_context(row, sample, stage, fieldname, input_file)}; "
                    "parallel braced field length mismatch: "
                    f"reference_field={reference_field}, "
                    f"expected={expected_length}, actual={actual_length}"
                )
    return reference_values


def _parse_braced_ints(row, fieldname, sample, stage, input_file):
    values = _validate_parallel_braced_fields(
        row,
        fieldname,
        sample,
        stage,
        input_file,
    )
    parsed_values = []
    for value in values:
        try:
            parsed_values.append(int(value))
        except (TypeError, ValueError) as exc:
            raise ValueError(
                f"{_row_error_context(row, sample, stage, fieldname, input_file)}; "
                f"invalid integer value={value!r}"
            ) from exc
    return parsed_values


def filter_levenshtein_row(
    row,
    max_distance,
    sample="<unknown>",
    stage="[2/8]",
    input_file="<unknown>",
):
    distances = _parse_braced_ints(
        row,
        "levenshtein_distance",
        sample,
        stage,
        input_file,
    )
    filtered_indices = [
        index for index, distance in enumerate(distances)
        if distance <= max_distance
    ]
    if not filtered_indices:
        return None

    new_row = {}
    for fieldname, value in row.items():
        if fieldname == "levenshtein_distance":
            new_row[fieldname] = join_braced(
                distances[index] for index in filtered_indices
            )
        elif isinstance(value, str) and "{" in value and "}" in value:
            entries = split_braced(value)
            new_row[fieldname] = join_braced(
                entries[index] for index in filtered_indices
            )
        else:
            new_row[fieldname] = value

    new_row["match_count"] = len(filtered_indices)
    new_row["UMI_count"] = len(set(split_braced(new_row["UMI"])))
    new_row["aaSeqCDR3_count"] = len(
        set(split_braced(new_row["aaSeqCDR3"]))
    )
    new_row["cloneid_count"] = len(set(split_braced(new_row["cloneId"])))
    return new_row


def add_umi_support_row(
    row,
    sample="<unknown>",
    stage="[3/8]",
    input_file="<unknown>",
):
    _validate_parallel_braced_fields(
        row,
        "UMI",
        sample,
        stage,
        input_file,
    )
    umis = re.findall(r"{([^}]+)}", row["UMI"])
    umi_counts = Counter(umis)
    new_row = dict(row)
    new_row["UMI_counts_private"] = join_braced(
        umi_counts[umi] for umi in umis
    )
    return new_row


def filter_umi_support_row(
    row,
    threshold,
    sample="<unknown>",
    stage="[4/8]",
    input_file="<unknown>",
):
    umi_counts_private = _parse_braced_ints(
        row,
        "UMI_counts_private",
        sample,
        stage,
        input_file,
    )
    filtered_indices = [
        index for index, count in enumerate(umi_counts_private)
        if count >= threshold
    ]
    if not filtered_indices:
        return None

    new_row = {}
    for fieldname, value in row.items():
        if fieldname == "UMI_counts_private":
            new_row[fieldname] = join_braced(
                umi_counts_private[index] for index in filtered_indices
            )
        elif isinstance(value, str) and "{" in value and "}" in value:
            entries = split_braced(value)
            new_row[fieldname] = join_braced(
                entries[index] for index in filtered_indices
            )
        else:
            new_row[fieldname] = value

    new_row["match_count"] = len(filtered_indices)
    new_row["UMI_count"] = len(set(split_braced(new_row["UMI"])))
    new_row["aaSeqCDR3_count"] = len(
        set(split_braced(new_row["aaSeqCDR3"]))
    )
    new_row["cloneid_count"] = len(set(split_braced(new_row["cloneId"])))
    return new_row


def _iter_rows_with_spatial_bc(input_file):
    rows = iter(iter_csv_rows(input_file))
    try:
        while True:
            try:
                row = next(rows)
            except StopIteration:
                break
            updated_row = dict(row)
            updated_row["spatial_bc"] = extract_unique_barcode(
                row["spatial_barcode_ori"]
            )
            del row
            yield updated_row
            del updated_row
    finally:
        close_rows = getattr(rows, "close", None)
        if close_rows is not None:
            close_rows()


def _iter_levenshtein_rows(input_file, max_distance, sample):
    rows = iter(iter_csv_rows(input_file))
    try:
        while True:
            try:
                row = next(rows)
            except StopIteration:
                break
            output_row = filter_levenshtein_row(
                row,
                max_distance,
                sample=sample,
                input_file=input_file,
            )
            del row
            yield output_row
            del output_row
    finally:
        close_rows = getattr(rows, "close", None)
        if close_rows is not None:
            close_rows()


def _iter_umi_support_rows(input_file, sample):
    rows = iter(iter_csv_rows(input_file))
    try:
        while True:
            try:
                row = next(rows)
            except StopIteration:
                break
            output_row = add_umi_support_row(
                row,
                sample=sample,
                input_file=input_file,
            )
            del row
            yield output_row
            del output_row
    finally:
        close_rows = getattr(rows, "close", None)
        if close_rows is not None:
            close_rows()


def _iter_filtered_umi_rows(input_file, threshold, sample):
    rows = iter(iter_csv_rows(input_file))
    try:
        while True:
            try:
                row = next(rows)
            except StopIteration:
                break
            output_row = filter_umi_support_row(
                row,
                threshold,
                sample=sample,
                input_file=input_file,
            )
            del row
            yield output_row
            del output_row
    finally:
        close_rows = getattr(rows, "close", None)
        if close_rows is not None:
            close_rows()


def _write_lv_summary(input_file, filtered_files):
    count_fields = [
        "UMI_count",
        "match_count",
        "aaSeqCDR3_count",
        "cloneid_count",
    ]
    output_fields = ["spatial_bc"] + [
        f"{fieldname}_Lv1" for fieldname in count_fields
    ]
    other_distances = [
        distance for distance in max_levenshtein_distances if distance != 1
    ]
    for distance in other_distances:
        output_fields.extend(
            f"{fieldname}_Lv{distance}" for fieldname in count_fields
        )

    other_iterators = {
        distance: iter(iter_csv_rows(filtered_files[distance]))
        for distance in other_distances
    }
    current_rows = {
        distance: next(iterator, None)
        for distance, iterator in other_iterators.items()
    }

    def summary_rows():
        reference_rows = iter(iter_csv_rows(filtered_files[1]))
        try:
            while True:
                try:
                    reference_row = next(reference_rows)
                except StopIteration:
                    break
                barcode = reference_row["spatial_bc"]
                output_row = {"spatial_bc": barcode}
                for fieldname in count_fields:
                    output_row[f"{fieldname}_Lv1"] = reference_row[fieldname]

                for distance in other_distances:
                    current_row = current_rows[distance]
                    iterator = other_iterators[distance]
                    while (
                        current_row is not None
                        and current_row["spatial_bc"] < barcode
                    ):
                        current_rows[distance] = None
                        del current_row
                        current_row = next(iterator, None)
                    current_rows[distance] = current_row
                    for fieldname in count_fields:
                        output_row[f"{fieldname}_Lv{distance}"] = (
                            current_row[fieldname]
                            if current_row is not None
                            and current_row["spatial_bc"] == barcode
                            else ""
                        )
                del reference_row
                yield output_row
                del output_row
        finally:
            close_reference_rows = getattr(reference_rows, "close", None)
            if close_reference_rows is not None:
                close_reference_rows()
            for iterator in other_iterators.values():
                close_iterator = getattr(iterator, "close", None)
                if close_iterator is not None:
                    close_iterator()

    output_file = f"{_pipeline_base(input_file)}_TRB_lv_UMI_reads_counts.csv"
    _write_csv_rows_atomic(output_file, output_fields, summary_rows())


def add_levenshtein_outputs(input_files):
    print("[2/8] Applying Levenshtein-distance filter and writing Lv summaries")
    for input_file in input_files:
        sample = _sample_label(input_file)
        input_fields = _csv_fieldnames(input_file)
        updated_fields = _append_fieldnames(input_fields, ["spatial_bc"])
        _write_csv_rows_atomic(
            input_file,
            updated_fields,
            _iter_rows_with_spatial_bc(input_file),
        )
        print(f"Added 'spatial_bc' column and saved updated file: {input_file}")

        filtered_fields = _append_fieldnames(
            updated_fields,
            ["UMI_count", "aaSeqCDR3_count", "cloneid_count"],
        )
        filtered_files = {}
        for max_distance in max_levenshtein_distances:
            output_file = (
                f"{_pipeline_base(input_file)}"
                f"_metadata_with_levenshtein_region_{max_distance}.csv"
            )
            filtered_files[max_distance] = output_file
            print(f"Output: {output_file}")
            _write_csv_rows_atomic(
                output_file,
                filtered_fields,
                _iter_levenshtein_rows(input_file, max_distance, sample),
                write_header_if_empty=False,
            )

        _write_lv_summary(input_file, filtered_files)
        print(
            "Merged data saved to "
            f"{_pipeline_base(input_file)}_TRB_lv_UMI_reads_counts.csv"
        )


def add_umi_private_counts(input_files):
    print("[3/8] Counting UMI support within each spatial barcode")
    for umi_threshold in umi_thresholds:
        for max_distance in max_levenshtein_distances:
            for input_file in input_files:
                sample = _sample_label(input_file)
                file_path = (
                    f"{_pipeline_base(input_file)}"
                    f"_metadata_with_levenshtein_region_{max_distance}.csv"
                )
                output_file = (
                    f"{_pipeline_base(input_file)}"
                    "_metadata_with_levenshtein_region_"
                    f"umi{umi_threshold}_{max_distance}.csv"
                )
                output_fields = _append_fieldnames(
                    _csv_fieldnames(file_path),
                    ["UMI_counts_private"],
                )
                print(f"Output: {output_file}")
                _write_csv_rows_atomic(
                    output_file,
                    output_fields,
                    _iter_umi_support_rows(file_path, sample),
                )


def make_filtered_umi_files(input_files):
    print("[4/8] Filtering entries by UMI support threshold")
    for umi_threshold in umi_thresholds:
        for max_distance in max_levenshtein_distances:
            for input_file in input_files:
                sample = _sample_label(input_file)
                input_file_path = (
                    f"{_pipeline_base(input_file)}"
                    "_metadata_with_levenshtein_region_"
                    f"umi{umi_threshold}_{max_distance}.csv"
                )
                output_file = (
                    f"{_pipeline_base(input_file)}"
                    "_filtered_metadata_with_levenshtein_region_"
                    f"umi{umi_threshold}_{max_distance}.csv"
                )
                print(f"Output: {output_file}")
                _write_csv_rows_atomic(
                    output_file,
                    _csv_fieldnames(input_file_path),
                    _iter_filtered_umi_rows(
                        input_file_path,
                        umi_threshold,
                        sample,
                    ),
                    write_header_if_empty=False,
                )


def build_clone_cdr3_counts(filtered_path, connection):
    first_order = 0
    with connection:
        connection.execute("DROP TABLE IF EXISTS clone_cdr3_counts")
        connection.execute(
            """
            CREATE TABLE clone_cdr3_counts (
                clone_id TEXT NOT NULL,
                cdr3 TEXT NOT NULL,
                support INTEGER NOT NULL,
                first_order INTEGER NOT NULL,
                PRIMARY KEY (clone_id, cdr3)
            )
            """
        )
        for row in iter_csv_rows(filtered_path):
            clone_ids = split_braced(row["cloneId"])
            cdr3_values = split_braced(row["aaSeqCDR3"])
            if len(clone_ids) != len(cdr3_values):
                raise ValueError(
                    "parallel braced field length mismatch: "
                    f"cloneId={len(clone_ids)}, aaSeqCDR3={len(cdr3_values)}"
                )

            for clone_id, cdr3 in zip(clone_ids, cdr3_values):
                connection.execute(
                    """
                    INSERT INTO clone_cdr3_counts (
                        clone_id, cdr3, support, first_order
                    )
                    VALUES (?, ?, 1, ?)
                    ON CONFLICT (clone_id, cdr3) DO UPDATE SET
                        support = support + 1,
                        first_order = MIN(first_order, excluded.first_order)
                    """,
                    (clone_id, cdr3, first_order),
                )
                first_order += 1


def select_clone_cdr3_representatives(connection):
    representatives = {}
    cursor = connection.execute(
        """
        SELECT clone_id, cdr3
        FROM clone_cdr3_counts
        ORDER BY support DESC, first_order ASC
        """
    )
    try:
        for clone_id, cdr3 in cursor:
            representatives.setdefault(clone_id, cdr3)
    finally:
        cursor.close()
    return representatives


def resolve_barcode_clones(row, representatives):
    expanded_fields = {
        fieldname: split_braced(value)
        for fieldname, value in row.items()
        if isinstance(value, str) and "{" in value and "}" in value
    }
    required_fields = {
        "spatial_barcode_ori",
        "UMI",
        "aaSeqCDR3",
        "cloneId",
        "UMI_counts_private",
    }
    missing_fields = required_fields - expanded_fields.keys()
    if missing_fields:
        raise KeyError(f"missing expanded fields: {sorted(missing_fields)}")

    entry_count = len(expanded_fields["UMI"])
    for fieldname, values in expanded_fields.items():
        if len(values) != entry_count:
            raise ValueError(
                "parallel braced field length mismatch: "
                f"UMI={entry_count}, {fieldname}={len(values)}"
            )

    expanded = []
    for position in range(entry_count):
        try:
            support = int(expanded_fields["UMI_counts_private"][position])
        except (TypeError, ValueError):
            continue
        entry = {
            fieldname: values[position]
            for fieldname, values in expanded_fields.items()
        }
        entry["support"] = support
        expanded.append(entry)

    clone_support = {}
    for entry in expanded:
        key = (
            entry["spatial_barcode_ori"],
            entry["UMI"],
            entry["cloneId"],
        )
        clone_support[key] = clone_support.get(key, 0) + entry["support"]

    selected_clones = {}
    for key in sorted(clone_support):
        barcode_umi = key[:2]
        current = selected_clones.get(barcode_umi)
        if current is None or clone_support[key] > clone_support[current]:
            selected_clones[barcode_umi] = key

    filtered = [
        entry
        for entry in expanded
        if selected_clones[(entry["spatial_barcode_ori"], entry["UMI"])][2]
        == entry["cloneId"]
    ]
    if not filtered:
        return None

    output = {"spatial_barcode_ori": filtered[0]["spatial_barcode_ori"]}
    source_index = row["index"]
    output["index"] = join_braced(source_index for _entry in filtered)
    for fieldname in expanded_fields:
        if fieldname != "spatial_barcode_ori":
            output[fieldname] = join_braced(
                entry[fieldname] for entry in filtered
            )

    representative_cdr3 = [
        representatives[entry["cloneId"]] for entry in filtered
    ]
    output["aaSeqCDR3_most"] = join_braced(representative_cdr3)
    output["match_count"] = len(filtered)
    output["UMI_count"] = len({entry["UMI"] for entry in filtered})
    output["aaSeqCDR3_count"] = len(
        {entry["aaSeqCDR3"] for entry in filtered}
    )
    output["cloneid_count"] = len({entry["cloneId"] for entry in filtered})
    output["aaSeqCDR3_most_count"] = len(set(representative_cdr3))
    return output


def iter_resolved_clone_rows(filtered_path, representatives):
    rows = iter(iter_csv_rows(filtered_path))
    try:
        try:
            row = next(rows)
        except StopIteration:
            return

        expanded_fieldnames = [
            fieldname
            for fieldname, value in row.items()
            if isinstance(value, str) and "{" in value and "}" in value
        ]
        source_ordinal = 0
        while True:
            entry_count = len(split_braced(row.get("UMI", "")))
            resolution_row = {"index": source_ordinal}
            for fieldname in expanded_fieldnames:
                values = split_braced(row.get(fieldname, ""))
                if len(values) < entry_count:
                    values.extend("nan" for _ in range(entry_count - len(values)))
                resolution_row[fieldname] = join_braced(values[:entry_count])

            del row
            output_row = resolve_barcode_clones(
                resolution_row,
                representatives,
            )
            del resolution_row
            if output_row is not None:
                yield output_row
                del output_row

            source_ordinal += 1
            try:
                row = next(rows)
            except StopIteration:
                break
    finally:
        close_rows = getattr(rows, "close", None)
        if close_rows is not None:
            close_rows()


def write_resolved_clone_rows(
    filtered_path,
    output_path,
    representatives,
    connection,
):
    with connection:
        connection.execute("DROP TABLE IF EXISTS resolved_clone_rows")
        connection.execute(
            """
            CREATE TABLE resolved_clone_rows (
                barcode TEXT NOT NULL,
                source_order INTEGER NOT NULL,
                payload TEXT NOT NULL,
                PRIMARY KEY (barcode, source_order)
            )
            """
        )

        fieldnames = None
        source_order = 0
        resolved_rows = iter(
            iter_resolved_clone_rows(filtered_path, representatives)
        )
        try:
            while True:
                try:
                    row = next(resolved_rows)
                except StopIteration:
                    break
                if fieldnames is None:
                    fieldnames = list(row)
                payload = json.dumps(
                    [row.get(fieldname) for fieldname in fieldnames],
                    separators=(",", ":"),
                )
                connection.execute(
                    """
                    INSERT INTO resolved_clone_rows (
                        barcode, source_order, payload
                    )
                    VALUES (?, ?, ?)
                    """,
                    (row["spatial_barcode_ori"], source_order, payload),
                )
                source_order += 1
                del row
        finally:
            close_rows = getattr(resolved_rows, "close", None)
            if close_rows is not None:
                close_rows()

    if fieldnames is None:
        _write_csv_rows_atomic(
            output_path,
            [],
            iter(()),
            write_header_if_empty=False,
        )
        with connection:
            connection.execute("DROP TABLE resolved_clone_rows")
        return

    def output_rows():
        cursor = connection.execute(
            """
            SELECT payload
            FROM resolved_clone_rows
            ORDER BY barcode, source_order
            """
        )
        try:
            for (payload,) in cursor:
                yield dict(zip(fieldnames, json.loads(payload)))
        finally:
            cursor.close()

    _write_csv_rows_atomic(
        output_path,
        fieldnames,
        output_rows(),
        write_header_if_empty=False,
    )
    with connection:
        connection.execute("DROP TABLE resolved_clone_rows")


SUMMARY_COLUMNS = [
    "spatial_barcode_ori",
    "aaSeqCDR3_most_count",
    "aaSeqCDR3_count",
    "cloneid_count",
    "UMI_count",
    "match_count",
]
VDJC_COLUMNS = [
    "spatial_barcode_ori",
    "aaSeqCDR3_most",
    "aaSeqCDR3_most_count",
    "aaSeqCDR3_count",
    "bestVGene",
    "bestDGene",
    "bestJGene",
    "bestCGene",
    "allVAlignments",
    "allDAlignments",
    "allJAlignments",
    "allCAlignments",
    "UMI",
    "UMI_count",
    "nSeqCDR3",
]
GENE_UPDATE_COLUMNS = [
    "bestVGene",
    "bestDGene",
    "bestJGene",
    "bestCGene",
    "nSeqCDR3",
    "allVAlignments",
    "allDAlignments",
    "allJAlignments",
    "allCAlignments",
    "UMI",
]
MUTATION_ADDITIONS = [
    f"all{region}Alignments_{suffix}"
    for region in ("V", "D", "J", "C")
    for suffix in ("mutation_rate", "sequence_length", "mutation_count")
]
EXPANDED_MUTATION_COLUMNS = [
    "spatial_barcode_ori",
    "aaSeqCDR3_most",
    "nSeqCDR3",
    "bestVGene",
    "bestDGene",
    "bestJGene",
    "bestCGene",
    "UMI",
    "allVAlignments_mutation_rate",
    "allDAlignments_mutation_rate",
    "allJAlignments_mutation_rate",
    "allCAlignments_mutation_rate",
    "allVAlignments_sequence_length",
    "allDAlignments_sequence_length",
    "allJAlignments_sequence_length",
    "allCAlignments_sequence_length",
    "allVAlignments_mutation_count",
    "allDAlignments_mutation_count",
    "allJAlignments_mutation_count",
    "allCAlignments_mutation_count",
    "VDJ_mutation_rate",
    "VDJC_mutation_rate",
    "VJ_mutation_rate",
]
REGION_RATE_COLUMNS = (
    "VDJ_mutation_rate",
    "VDJC_mutation_rate",
    "VJ_mutation_rate",
)


class _ExpandedMutationWriter:
    def __init__(self, writer):
        self.writer = writer
        self.float_rate_fields = set()

    def writeheader(self):
        self.writer.writeheader()

    def writerow(self, row):
        for fieldname in REGION_RATE_COLUMNS:
            if isinstance(row[fieldname], float):
                self.float_rate_fields.add(fieldname)
        self.writer.writerow(row)


def _normalize_expanded_rate_zeros(path, float_rate_fields):
    if not float_rate_fields:
        return
    path = Path(path)
    descriptor, normalized_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".normalized.tmp",
    )
    os.close(descriptor)
    normalized_path = Path(normalized_name)
    with path.open(newline="", encoding="utf-8") as input_handle:
        reader = csv.DictReader(input_handle)
        with normalized_path.open("w", newline="", encoding="utf-8") as output_handle:
            writer = csv.DictWriter(
                output_handle,
                fieldnames=reader.fieldnames,
                lineterminator=os.linesep,
            )
            writer.writeheader()
            for row in reader:
                for fieldname in float_rate_fields:
                    if row[fieldname] == "0":
                        row[fieldname] = "0.0"
                writer.writerow(row)
    normalized_path.replace(path)


def _clone_output_path(input_file, umi_threshold, max_distance):
    return (
        f"{_pipeline_base(input_file)}"
        "_cloneid_metadata_with_levenshtein_region_"
        f"umi{umi_threshold}_lv{max_distance}.csv"
    )


def _iter_selected_rows(input_file, fieldnames):
    for row in iter_csv_rows(input_file):
        yield {fieldname: row[fieldname] for fieldname in fieldnames}


def make_filtered_summary_files(input_files):
    print("[6/8] Writing compact barcode-level count summary tables")
    compact_files = {input_file: [] for input_file in input_files}
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_distance in max_levenshtein_distances:
                source_file = _clone_output_path(
                    input_file, umi_threshold, max_distance
                )
                output_file = (
                    f"filtered_data_{_pipeline_base(input_file)}"
                    f"_umi{umi_threshold}_{max_distance}.csv"
                )
                _write_csv_rows_atomic(
                    output_file,
                    SUMMARY_COLUMNS,
                    _iter_selected_rows(source_file, SUMMARY_COLUMNS),
                )
                compact_files[input_file].append(
                    (umi_threshold, max_distance, output_file)
                )
                print(f"Filtered data saved to {output_file}")

    for input_file in input_files:
        merged_frame = None
        for umi_threshold, max_distance, file_path in compact_files[input_file]:
            frame = pd.read_csv(file_path).add_suffix(
                f"_Lv{max_distance}_UMI{umi_threshold}"
            )
            frame.rename(
                columns={
                    f"spatial_barcode_ori_Lv{max_distance}_UMI{umi_threshold}":
                    "spatial_barcode_ori"
                },
                inplace=True,
            )
            merged_frame = (
                frame
                if merged_frame is None
                else pd.merge(
                    merged_frame,
                    frame,
                    on="spatial_barcode_ori",
                    how="outer",
                )
            )
        output_file = f"merged_filtered_data_{_pipeline_base(input_file)}.csv"
        if merged_frame is None:
            merged_frame = pd.DataFrame()
        _write_frame_atomic(
            merged_frame,
            output_file,
            index=False,
            encoding="utf-8",
        )
        del merged_frame
        print(f"Final merged data saved to {output_file}")


def _write_matrix_combination(input_file, umi_threshold, max_distance):
    source_file = _clone_output_path(input_file, umi_threshold, max_distance)
    cdr3_counts = CounterState()
    umi_sets = CounterState()
    for row in iter_csv_rows(source_file):
        update_matrix_counters(row, cdr3_counts, umi_sets)

    cdr3_frame = pd.DataFrame.from_dict(cdr3_counts, orient="index").fillna(0)
    umi_counts = {
        barcode: {
            cdr3: len(umis)
            for cdr3, umis in barcode_values.items()
        }
        for barcode, barcode_values in umi_sets.items()
    }
    umi_frame = pd.DataFrame.from_dict(umi_counts, orient="index").fillna(0)
    cdr3_output = (
        f"{_pipeline_base(input_file)}_aaSeqCDR3_ori_count_matrix_"
        f"umi{umi_threshold}_levenshtein_{max_distance}.csv"
    )
    umi_output = (
        f"{_pipeline_base(input_file)}_aaSeqCDR3_UMI_count_matrix_"
        f"umi{umi_threshold}_levenshtein_{max_distance}.csv"
    )
    with atomic_output_group([cdr3_output, umi_output]) as temporary_paths:
        cdr3_frame.to_csv(
            temporary_paths[0], index=True, lineterminator=os.linesep
        )
        umi_frame.to_csv(
            temporary_paths[1], index=True, lineterminator=os.linesep
        )

    label = (
        f"{_pipeline_base(input_file)}_umi{umi_threshold}"
        f"_levenshtein_{max_distance}"
    )
    cdr3_summary = {
        "file_name": label,
        "aaSeqCDR3_count": cdr3_frame.shape[1],
        "spatial_barcode_count": cdr3_frame.shape[0],
    }
    umi_summary = {
        "file_name": label,
        "aaSeqCDR3_count": umi_frame.shape[1],
        "spatial_barcode_count": umi_frame.shape[0],
    }
    del cdr3_frame, umi_frame, umi_counts, cdr3_counts, umi_sets
    return cdr3_summary, umi_summary


def write_matrix_summaries(input_files):
    print("[7/8] Generating aaSeqCDR3 count and UMI count matrices")
    cdr3_summaries = []
    umi_summaries = []
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_distance in max_levenshtein_distances:
                cdr3_summary, umi_summary = _write_matrix_combination(
                    input_file, umi_threshold, max_distance
                )
                cdr3_summaries.append(cdr3_summary)
                umi_summaries.append(umi_summary)

    cdr3_summary_frame = pd.DataFrame(cdr3_summaries)
    _write_frame_atomic(
        cdr3_summary_frame,
        f"{_pipeline_base(input_files[-1])}_unique_aaSeqCDR3_count.csv",
        index=False,
        encoding="utf-8",
    )
    print(cdr3_summary_frame)
    umi_summary_frame = pd.DataFrame(umi_summaries)
    _write_frame_atomic(
        umi_summary_frame,
        f"{_pipeline_base(input_files[-1])}_unique_UMI_aaSeqCDR3_count.csv",
        index=False,
        encoding="utf-8",
    )
    print(umi_summary_frame)
    del cdr3_summary_frame, umi_summary_frame


def write_shannon(input_files):
    print("[8/8] Calculating Shannon diversity and VDJC/mutation summaries")
    all_shannon_diversities = pd.DataFrame()
    last_input_file = input_files[-1]
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            last_input_file = input_file
            for max_distance in max_levenshtein_distances:
                matrix_file = (
                    f"{_pipeline_base(input_file)}"
                    "_aaSeqCDR3_UMI_count_matrix_"
                    f"umi{umi_threshold}_levenshtein_{max_distance}.csv"
                )
                matrix_frame = pd.read_csv(matrix_file, index_col=0)
                shannon_values = matrix_frame.apply(
                    lambda row: calculate_shannon_from_counts(row.values),
                    axis=1,
                )
                column_name = (
                    f"Shannon_count_Lv{max_distance}_UMI{umi_threshold}"
                )
                if all_shannon_diversities.empty:
                    all_shannon_diversities = pd.DataFrame(
                        shannon_values, columns=[column_name]
                    )
                else:
                    all_shannon_diversities[column_name] = shannon_values
                del matrix_frame, shannon_values

    _write_frame_atomic(
        all_shannon_diversities,
        f"{_pipeline_base(last_input_file)}_merged_shannon_diversity.csv",
        index=True,
    )
    del all_shannon_diversities
    print("Shannon diversity index calculation and merging completed.")


def _update_genes_for_barcode(row, sample, input_file):
    _validate_parallel_braced_fields(
        row,
        "UMI",
        sample,
        "[8/8] VDJC",
        input_file,
    )
    umis = split_braced(row["UMI"])
    parsed = {
        fieldname: split_braced(row[fieldname])
        for fieldname in GENE_UPDATE_COLUMNS
        if fieldname != "UMI"
    }
    updated = {fieldname: [""] * len(umis) for fieldname in parsed}
    umi_indices = {}
    for index, umi in enumerate(umis):
        umi_indices.setdefault(umi, []).append(index)
    for indices in umi_indices.values():
        for fieldname, values in parsed.items():
            representative = Counter(values[index] for index in indices).most_common(1)[0][0]
            for index in indices:
                updated[fieldname][index] = representative

    output = dict(row)
    for fieldname, values in updated.items():
        output[fieldname] = join_braced(values)
    return output


def _write_vdjc_combination(input_file, umi_threshold, max_distance):
    source_file = _clone_output_path(input_file, umi_threshold, max_distance)
    sample = _sample_label(input_file)
    source_fields = _csv_fieldnames(source_file)
    mutation_fields = _append_fieldnames(source_fields, MUTATION_ADDITIONS)
    base = _pipeline_base(input_file)
    output_files = [
        f"{base}_VDJC_umi{umi_threshold}_lv{max_distance}.csv",
        f"{base}_VDJC_VJ_frequency_umi{umi_threshold}_lv{max_distance}.csv",
        f"{base}_with_mutation_rates_umi{umi_threshold}_lv{max_distance}.csv",
        f"{base}_expanded_with_region_rates_umi{umi_threshold}_lv{max_distance}.csv",
    ]

    gene_umi_sets = CounterState()
    vj_umi_sets = CounterState()
    with atomic_output_group(output_files) as temporary_paths:
        with contextlib.ExitStack() as stack:
            handles = [
                stack.enter_context(
                    path.open("w", newline="", encoding="utf-8")
                )
                for path in temporary_paths
            ]
            vdjc_writer = csv.DictWriter(
                handles[0], fieldnames=VDJC_COLUMNS, lineterminator=os.linesep
            )
            frequency_writer = csv.DictWriter(
                handles[1],
                fieldnames=["Gene", "Frequency", "Category"],
                lineterminator=os.linesep,
            )
            mutation_writer = csv.DictWriter(
                handles[2], fieldnames=mutation_fields, lineterminator=os.linesep
            )
            expanded_writer = _ExpandedMutationWriter(
                csv.DictWriter(
                    handles[3],
                    fieldnames=EXPANDED_MUTATION_COLUMNS,
                    lineterminator=os.linesep,
                )
            )
            for writer in (
                vdjc_writer,
                frequency_writer,
                mutation_writer,
                expanded_writer,
            ):
                writer.writeheader()

            for row in iter_csv_rows(source_file):
                updated_row = _update_genes_for_barcode(
                    row, sample, source_file
                )
                vdjc_writer.writerow(
                    {
                        fieldname: updated_row[fieldname]
                        for fieldname in VDJC_COLUMNS
                    }
                )
                update_gene_frequency_counters(
                    updated_row, gene_umi_sets, vj_umi_sets
                )
                write_mutation_rows(
                    updated_row,
                    mutation_writer,
                    expanded_writer,
                    sample=sample,
                    stage="[8/8] mutation",
                    input_file=source_file,
                )
                del row, updated_row

            for category, gene_sets in gene_umi_sets.items():
                for gene, umis in gene_sets.items():
                    frequency_writer.writerow(
                        {
                            "Gene": gene,
                            "Frequency": len(umis),
                            "Category": category,
                        }
                    )
            for pair, umis in vj_umi_sets.items():
                frequency_writer.writerow(
                    {
                        "Gene": pair,
                        "Frequency": len(umis),
                        "Category": "V-J Paired Gene",
                    }
                )
            float_rate_fields = set(expanded_writer.float_rate_fields)
        del handles, vdjc_writer, frequency_writer, mutation_writer, expanded_writer
        _normalize_expanded_rate_zeros(
            temporary_paths[3], float_rate_fields
        )
        del float_rate_fields
    del gene_umi_sets, vj_umi_sets
    print(
        f"VDJC and V-J paired gene frequency calculation completed for "
        f"{source_file}"
    )
    print(
        f"Processed {source_file}. Mutation rates saved: {output_files[2]}, "
        f"Expanded saved: {output_files[3]}"
    )


def write_vdjc_outputs(input_files):
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_distance in max_levenshtein_distances:
                _write_vdjc_combination(
                    input_file, umi_threshold, max_distance
                )


def run_summary_pipeline(input_files):
    for input_file in input_files:
        for umi_threshold in umi_thresholds:
            for max_distance in max_levenshtein_distances:
                source_file = _clone_output_path(
                    input_file, umi_threshold, max_distance
                )
                if not Path(source_file).is_file():
                    raise FileNotFoundError(
                        f"sample={_sample_label(input_file)}; stage=[6/8]-[8/8]; "
                        f"input_file={source_file}; prior-stage cloneId output is missing"
                    )
    make_filtered_summary_files(input_files)
    write_matrix_summaries(input_files)
    write_shannon(input_files)
    write_vdjc_outputs(input_files)


def find_8A_position_and_generate_windows(sequence):
    # Candidate barcodes start at base 6 and occupy 16 bp windows in the next 27 bp.
    start_index = 6
    end_index = min(6 + 27, len(sequence))
    return [sequence[i:i + 16] for i in range(start_index, end_index - 15)]


def extract_umi(target_sequence_part, original_sequence):
    # UMI is the 12 bp immediately upstream of the observed barcode.
    start_index = target_sequence_part.find(original_sequence)
    if start_index == -1:
        return ""
    umi = target_sequence_part[start_index - 12:start_index] if start_index >= 12 else target_sequence_part[:start_index]
    return umi.rjust(12, "A")


def build_substitution_index(whitelist_values):
    """Validate 16 bp whitelist barcodes and index every unique substitution."""
    values = list(whitelist_values)
    if not values:
        raise ValueError("spatial barcode whitelist is empty")
    if len(values) != len(set(values)):
        raise ValueError("spatial barcode whitelist contains duplicate values")

    whitelist = set(values)
    for barcode in whitelist:
        if not isinstance(barcode, str) or len(barcode) != 16:
            raise ValueError(f"invalid spatial barcode length: {barcode!r}")
        invalid_bases = set(barcode) - set(DNA_BASES)
        if invalid_bases:
            raise ValueError(
                f"invalid spatial barcode bases in {barcode!r}: "
                f"{''.join(sorted(invalid_bases))}"
            )

    substitution_index = {}
    conflicts = {}
    for barcode in sorted(whitelist):
        for position, original_base in enumerate(barcode):
            for replacement in DNA_BASES:
                if replacement == original_base:
                    continue
                neighbor = (
                    barcode[:position]
                    + replacement
                    + barcode[position + 1:]
                )
                if neighbor in whitelist:
                    continue
                candidate = (barcode, position)
                existing = substitution_index.get(neighbor)
                if existing is None:
                    substitution_index[neighbor] = candidate
                elif existing != candidate:
                    conflicts.setdefault(neighbor, {existing}).add(candidate)

    if conflicts:
        examples = []
        for neighbor in sorted(conflicts)[:3]:
            candidates = sorted(barcode for barcode, _ in conflicts[neighbor])
            examples.append(f"{neighbor}=>{','.join(candidates)}")
        raise ValueError(
            "ambiguous single-substitution barcode neighbors detected: "
            f"count={len(conflicts)}; examples={'; '.join(examples)}"
        )

    return whitelist, substitution_index


def process_sequence(row, whitelist, substitution_index):
    """Match one MiXCR alignment row to a spatial barcode and extract UMI."""
    try:
        target_sequence = row["targetSequences"].split(",")[0]
        target_qualities = row["targetQualities"][:len(target_sequence)]

        windows = find_8A_position_and_generate_windows(target_sequence)
        if not windows:
            return None

        quality_windows = [
            target_qualities[6 + i:6 + i + 16]
            for i in range(len(windows))
        ]

        matched_qualities = []
        operation_types = []

        for i, window in enumerate(windows):
            if window in whitelist:
                matched_qualities.append(quality_windows[i])
                operation_types.append("match")
                umi = extract_umi(target_sequence, window)
                updated_row = dict(row)
                updated_row["original_sequence"] = window
                updated_row["corrected_sequence"] = window
                updated_row["matched_spatial_barcode"] = window
                updated_row["match_count"] = 1
                updated_row["matched_qualities"] = ",".join(matched_qualities)
                updated_row["levenshtein_distance"] = 0
                updated_row["operation_type"] = ",".join(operation_types)
                updated_row["UMI"] = umi
                return updated_row

        for i, window in enumerate(windows):
            indexed_match = substitution_index.get(window)
            if indexed_match is None:
                continue
            whitelist_seq, mismatch_position = indexed_match
            if (ord(quality_windows[i][mismatch_position]) - 33) >= 30:
                continue

            best_match = {
                "original_sequence": window,
                "corrected_sequence": whitelist_seq,
                "matched_spatial_barcode": whitelist_seq,
                "match_count": 1,
                "matched_qualities": quality_windows[i],
                "levenshtein_distance": 1,
                "operation_type": "mismatch",
            }
            umi = extract_umi(target_sequence, window)
            updated_row = dict(row)
            updated_row.update(best_match)
            updated_row["UMI"] = umi
            return updated_row
    except Exception as exc:
        print(f"Error while processing row: {exc}")
    return None


def process_alignments_to_store(
    input_file, metadata, store, chunksize, output_file=None
):
    """Match alignment chunks and commit each chunk's matches to SQLite."""
    if "spatial_barcode_mixcr" not in metadata.columns:
        raise ValueError("missing spatial_barcode_mixcr column in metadata")
    whitelist, substitution_index = build_substitution_index(
        metadata["spatial_barcode_mixcr"].tolist()
    )
    print(
        "Built spatial barcode substitution index: "
        f"whitelist={len(whitelist)}, neighbors={len(substitution_index)}"
    )
    if output_file is not None:
        print(f"Input:  {input_file}")
        print(f"Output: {output_file}")
    scanned_rows = 0
    matched_rows = 0
    matched_barcodes = set()

    for chunk in pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize):
        scanned_rows += len(chunk)
        updated_records = []
        for row in chunk.replace("", float("nan")).to_dict("records"):
            result = process_sequence(row, whitelist, substitution_index)
            if result:
                updated_records.append(result)

        matched_rows += len(updated_records)
        matched_barcodes.update(
            record["matched_spatial_barcode"] for record in updated_records
        )

        if not updated_records:
            print(f"Scanned {scanned_rows} rows; no matches in current chunk")
            continue

        updated_records = (
            pd.DataFrame(updated_records)
            .merge(
                metadata,
                left_on="matched_spatial_barcode",
                right_on="spatial_barcode_mixcr",
                how="left",
            )
            .to_dict("records")
        )
        if list(updated_records[0]) != store.columns:
            raise ValueError("store columns do not match stage 1 output columns")
        store.add_records(updated_records, {})
        print(
            f"Scanned {scanned_rows} rows; matched {matched_rows} rows; "
            f"barcode groups: {len(matched_barcodes)}"
        )
        del updated_records


def write_stage1_aggregated(store, output_file, output_columns):
    """Stream one aggregated CSV row per barcode in deterministic order."""
    with atomic_output(output_file) as temporary_path:
        wrote_header = False
        groups = iter(store.iter_barcode_groups())
        try:
            while True:
                try:
                    barcode, match_count, values_by_column = next(groups)
                except StopIteration:
                    break

                row = {}
                for column in output_columns:
                    if column == "match_count":
                        row[column] = match_count
                    else:
                        row[column] = (
                            "{"
                            + "}{".join(
                                map(str, values_by_column.get(column, []))
                            )
                            + "}"
                        )
                pd.DataFrame([row], columns=output_columns).to_csv(
                    temporary_path,
                    mode="a" if wrote_header else "w",
                    header=not wrote_header,
                    index=False,
                )
                wrote_header = True

                del row
                del values_by_column
                del match_count
                del barcode
        finally:
            close_groups = getattr(groups, "close", None)
            if close_groups is not None:
                close_groups()

        if not wrote_header:
            pd.DataFrame().to_csv(temporary_path, index=False)
        return wrote_header


class DiskMatchStore:
    """Store matched reads on disk and stream them in deterministic groups."""

    def __init__(self, path, columns):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.columns = list(columns)
        self.value_columns = [column for column in self.columns if column != "match_count"]
        self.connection = sqlite3.connect(self.path)
        self.connection.execute("PRAGMA journal_mode=OFF")
        self.connection.execute("PRAGMA synchronous=OFF")
        self.connection.execute("PRAGMA temp_store=FILE")
        self.connection.execute(
            """
            CREATE TABLE matched_reads (
                read_order INTEGER PRIMARY KEY,
                barcode TEXT NOT NULL,
                match_count INTEGER NOT NULL,
                payload TEXT NOT NULL
            )
            """
        )
        self.connection.execute(
            """
            CREATE INDEX matched_reads_barcode_order
            ON matched_reads(barcode, read_order)
            """
        )
        self.connection.commit()

    def add_records(self, records, metadata_by_barcode):
        """Insert one input chunk in a single transaction."""
        rows = []
        for record in records:
            barcode = record["matched_spatial_barcode"]
            values = {**metadata_by_barcode.get(barcode, {}), **record}
            payload = [str(values.get(column, "")) for column in self.value_columns]
            rows.append(
                (
                    barcode,
                    int(record.get("match_count", 1)),
                    json.dumps(payload, separators=(",", ":")),
                )
            )

        with self.connection:
            self.connection.executemany(
                """
                INSERT INTO matched_reads (barcode, match_count, payload)
                VALUES (?, ?, ?)
                """,
                rows,
            )

    def iter_barcode_groups(self):
        """Yield barcode groups ordered by barcode and original read order."""
        cursor = self.connection.execute(
            """
            SELECT barcode, match_count, payload
            FROM matched_reads
            ORDER BY barcode, read_order
            """
        )
        barcode = None
        match_count = 0
        values_by_column = None

        try:
            for current_barcode, current_match_count, payload in cursor:
                if barcode is not None and current_barcode != barcode:
                    yield barcode, match_count, values_by_column
                    match_count = 0
                    values_by_column = None

                if values_by_column is None:
                    barcode = current_barcode
                    values_by_column = {column: [] for column in self.value_columns}

                match_count += current_match_count
                for column, value in zip(self.value_columns, json.loads(payload)):
                    values_by_column[column].append(value)

            if barcode is not None:
                yield barcode, match_count, values_by_column
        finally:
            cursor.close()

    def close(self):
        self.connection.close()


@contextlib.contextmanager
def atomic_output(final_path):
    """Expose a temporary path and atomically publish it after success."""
    final_path = Path(final_path)
    final_path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=final_path.parent,
        prefix=f".{final_path.name}.",
        suffix=".tmp",
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        yield temporary_path
    except BaseException:
        raise
    else:
        temporary_path.replace(final_path)


@contextlib.contextmanager
def atomic_output_group(final_paths):
    """Publish a related output set together, restoring old files on failure."""
    final_paths = [Path(path) for path in final_paths]
    if len(set(final_paths)) != len(final_paths):
        raise ValueError("atomic output paths must be unique")

    temporary_paths = []
    for final_path in final_paths:
        final_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            dir=final_path.parent,
            prefix=f".{final_path.name}.",
            suffix=".tmp",
        )
        os.close(descriptor)
        temporary_paths.append(Path(temporary_name))

    try:
        yield temporary_paths
    except BaseException:
        raise
    else:
        backups = {}
        published = []
        try:
            for final_path in final_paths:
                if final_path.exists():
                    descriptor, backup_name = tempfile.mkstemp(
                        dir=final_path.parent,
                        prefix=f".{final_path.name}.",
                        suffix=".bak",
                    )
                    os.close(descriptor)
                    backup_path = Path(backup_name)
                    backup_path.unlink()
                    final_path.replace(backup_path)
                    backups[final_path] = backup_path
            for temporary_path, final_path in zip(
                temporary_paths, final_paths
            ):
                temporary_path.replace(final_path)
                published.append(final_path)
        except BaseException:
            for final_path in published:
                if final_path.exists():
                    final_path.unlink()
            for final_path, backup_path in backups.items():
                if backup_path.exists():
                    backup_path.replace(final_path)
            raise
        else:
            for backup_path in backups.values():
                if backup_path.exists():
                    backup_path.unlink()


def run_diskbacked_pipeline(sample_name, base_path, chunksize):
    """Run all eight stages while retaining failed spool state for diagnosis."""
    base_path = Path(base_path)
    metadata_path = base_path / f"{sample_name}_meta_data_with_mixcr.csv"
    alignment_path = base_path / f"{sample_name}-filtered_alignments_poly8A.tsv"
    output_name = f"{sample_name}_updated_alignments_umi_poly8A3.csv"
    output_path = base_path / output_name
    sqlite_path = base_path / f".{sample_name}.stage6.sqlite3"
    recovery_artifacts = {sqlite_path}
    recovery_artifacts.update(
        candidate
        for candidate in base_path.iterdir()
        if (
            candidate.name.startswith(".")
            and sample_name in candidate.name
            and candidate.name.endswith((".tmp", ".bak"))
        )
    )
    retained_state = sorted(
        path for path in recovery_artifacts if path.exists()
    )
    if retained_state:
        raise RuntimeError(
            f"sample={sample_name} stage=preflight refusing to overwrite "
            "retained stage 6 recovery state: "
            + ", ".join(str(path) for path in retained_state)
        )

    metadata = pd.read_csv(metadata_path)
    alignment_fields = list(
        pd.read_csv(alignment_path, sep="\t", dtype=str, nrows=0).columns
    )
    generated_fields = [
        "original_sequence",
        "corrected_sequence",
        "matched_spatial_barcode",
        "match_count",
        "matched_qualities",
        "levenshtein_distance",
        "operation_type",
        "UMI",
    ]
    output_fields = list(
        pd.DataFrame(columns=[*alignment_fields, *generated_fields])
        .merge(
            metadata.iloc[0:0],
            left_on="matched_spatial_barcode",
            right_on="spatial_barcode_mixcr",
            how="left",
        )
        .columns
    )

    store = DiskMatchStore(sqlite_path, output_fields)
    try:
        print(
            "[1/8] Matching MiXCR alignments to spatial barcode whitelist "
            "and extracting UMI"
        )
        process_alignments_to_store(
            alignment_path,
            metadata,
            store,
            chunksize=chunksize,
            output_file=output_path,
        )
        has_matches = write_stage1_aggregated(
            store,
            output_path,
            output_fields,
        )
        if has_matches:
            print(f"Processing complete; results saved to '{output_path}'")
        else:
            print(f"No matches found; saved empty file: {output_path}")

        input_files = [output_name]
        add_levenshtein_outputs(input_files)
        add_umi_private_counts(input_files)
        make_filtered_umi_files(input_files)
        print("[5/8] Resolving cloneId per spatial barcode and UMI")
        for umi_threshold in umi_thresholds:
            for max_distance in max_levenshtein_distances:
                filtered_path = (
                    f"{_pipeline_base(output_name)}"
                    "_filtered_metadata_with_levenshtein_region_"
                    f"umi{umi_threshold}_{max_distance}.csv"
                )
                clone_path = _clone_output_path(
                    output_name, umi_threshold, max_distance
                )
                build_clone_cdr3_counts(filtered_path, store.connection)
                representatives = select_clone_cdr3_representatives(
                    store.connection
                )
                write_resolved_clone_rows(
                    filtered_path,
                    clone_path,
                    representatives,
                    store.connection,
                )
                print(f"Processed data saved to {clone_path}")
        run_summary_pipeline(input_files)
    except BaseException:
        store.close()
        raise
    else:
        store.close()
        sqlite_path.unlink()


if int(os.environ["MAX_WORKERS"]) != 1:
    raise ValueError("06_spatial_barcode_umi_pipeline_diskbacked.sh requires MAX_WORKERS=1")

sample_name = os.environ.get("SAMPLE_NAME")
if sample_name:
    base_path = Path(os.environ.get("MIXCR_DIR", "."))
    metadata_path = base_path / f"{sample_name}_meta_data_with_mixcr.csv"
    alignment_path = base_path / f"{sample_name}-filtered_alignments_poly8A.tsv"
    if os.environ.get("STAGE6_SUMMARY_ONLY") == "1":
        run_summary_pipeline(
            [f"{sample_name}_updated_alignments_umi_poly8A3.csv"]
        )
    else:
        missing_paths = [
            path
            for path in (metadata_path, alignment_path)
            if not path.is_file()
        ]
        if missing_paths:
            raise FileNotFoundError(
                f"sample={sample_name} stage=entry missing required input(s): "
                + ", ".join(str(path) for path in missing_paths)
                + "; set STAGE6_SUMMARY_ONLY=1 only for an explicit "
                "summary-only operation"
            )
        run_diskbacked_pipeline(
            sample_name,
            base_path,
            int(os.environ.get("CHUNKSIZE", "50000")),
        )
PY
    ); then
        processed_count=$((processed_count + 1))
    else
        sample_status=$?
        failed_count=$((failed_count + 1))
        echo "Stage 6 sample failed: sample=$sample_name exit=$sample_status" >&2
    fi
done

echo "Stage 6 sample summary: discovered=$discovered_count processed=$processed_count failed=$failed_count"
if [ "$failed_count" -ne 0 ] || [ "$processed_count" -eq 0 ]; then
    exit 1
fi
