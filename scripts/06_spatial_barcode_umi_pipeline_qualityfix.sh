#!/usr/bin/env bash
set -euo pipefail

# Low-memory rewrite of the original IGH/TRB pipeline.
# Calculation rules, thresholds, and output names are kept aligned with the
# original script; the main change is that the huge alignment TSV is processed
# in chunks instead of being loaded and queued in memory all at once.

BASE_DIR="${BASE_DIR:-results}"
SAMPLE_GLOB="${SAMPLE_GLOB:-*-TRB}"
PREPROCESS_DIR="${PREPROCESS_DIR:-preprocessing}"
MIXCR_OUTDIR="${MIXCR_OUTDIR:-mixcr}"
CHUNKSIZE="${CHUNKSIZE:-50000}"
MAX_WORKERS="${MAX_WORKERS:-1}"

# Process filtered MiXCR alignments sample by sample. CHUNKSIZE controls how
# many alignment rows are held in memory at once; MAX_WORKERS can be increased
# cautiously for faster barcode matching on larger machines.
find "$BASE_DIR" -maxdepth 1 -type d -name "$SAMPLE_GLOB" | sort | while read -r sample_dir; do
    sample_name=$(basename "$sample_dir")
    mixcr_dir="${sample_dir}/${PREPROCESS_DIR}/${MIXCR_OUTDIR}"

    echo "========================================"
    echo "Processing ${sample_name}"
    echo "Dir: ${mixcr_dir}"

    if [ ! -d "$mixcr_dir" ]; then
        echo "Missing directory, skipping: $mixcr_dir"
        continue
    fi

    if [ ! -f "${mixcr_dir}/${sample_name}_meta_data_with_mixcr.csv" ]; then
        echo "Missing metadata, skipping: ${mixcr_dir}/${sample_name}_meta_data_with_mixcr.csv"
        continue
    fi

    if [ ! -f "${mixcr_dir}/${sample_name}-filtered_alignments_poly8A.tsv" ]; then
        echo "Missing filtered alignments, skipping: ${mixcr_dir}/${sample_name}-filtered_alignments_poly8A.tsv"
        continue
    fi

    cd "$mixcr_dir"

    SAMPLE_NAME="$sample_name" \
    MIXCR_DIR="$mixcr_dir" \
    CHUNKSIZE="$CHUNKSIZE" \
    MAX_WORKERS="$MAX_WORKERS" \
    python3 - <<'PY'
"""Spatial barcode, UMI, clonotype, and VDJC summary workflow.

Inputs in MIXCR_DIR:
  - {sample_name}-filtered_alignments_poly8A.tsv:
    filtered MiXCR alignments from scripts/05_filter_alignments.sh.
  - {sample_name}_meta_data_with_mixcr.csv:
    spatial metadata containing the whitelist column spatial_barcode_mixcr.

Main outputs:
  - {sample_name}_updated_alignments_umi_poly8A3.csv
  - *_metadata_with_levenshtein_region_*.csv
  - *_filtered_metadata_with_levenshtein_region_umi*.csv
  - *_cloneid_metadata_with_levenshtein_region_umi*_lv*.csv
  - *_aaSeqCDR3_*_matrix_*.csv
  - *_VDJC_*.csv, *_VDJC_VJ_frequency_*.csv
  - *_with_mutation_rates_*.csv
  - *_expanded_with_region_rates_*.csv
"""
import math
import os
import re
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from rapidfuzz.distance import Levenshtein


sample_name = os.environ["SAMPLE_NAME"]
base_path = os.environ["MIXCR_DIR"]
chunksize = int(os.environ.get("CHUNKSIZE", "50000"))
max_workers = max(1, int(os.environ.get("MAX_WORKERS", "1")))

meta_data_files = {
    f"{sample_name}_meta_data_with_mixcr.csv": [sample_name],
}
max_levenshtein_distances = [1]
umi_thresholds = [3]


# Curly-brace strings store all read-level values belonging to one spatial
# barcode in one CSV cell, e.g. {UMI1}{UMI2}{UMI3}. These helpers preserve the
# original representation used by downstream steps.
def split_braced(value):
    if not isinstance(value, str) or value == "":
        return []
    return value.strip("{}").split("}{")


def join_braced(values):
    return "{" + "}{".join(map(str, values)) + "}"


def is_missing(value):
    return value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value)


def find_8A_position_and_generate_windows(sequence):
    # After polyA trimming, candidate spatial barcodes are searched in a fixed
    # window: start at base 6 and scan the next 27 bp with 16 bp windows.
    start_index = 6
    end_index = min(6 + 27, len(sequence))
    return [sequence[i:i + 16] for i in range(start_index, end_index - 15)]


def extract_umi(target_sequence_part, original_sequence):
    # UMI is defined as the 12 bp immediately upstream of the matched barcode.
    # If fewer than 12 bp are available, left-pad with A to preserve length.
    start_index = target_sequence_part.find(original_sequence)
    if start_index == -1:
        return ""
    umi = target_sequence_part[start_index - 12:start_index] if start_index >= 12 else target_sequence_part[:start_index]
    return umi.rjust(12, "A")


def process_sequence(row, whitelist):
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
            # First pass: exact whitelist match. This is the highest-confidence
            # case and keeps Levenshtein distance at 0.
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

        best_match = None
        # Second pass: allow edit distance 1. Mismatch corrections are accepted
        # only when the differing base quality is below Q30.
        for lev_dist in range(1, 2):
            for i, window in enumerate(windows):
                for whitelist_seq in whitelist:
                    if Levenshtein.distance(window, whitelist_seq) != lev_dist:
                        continue

                    ops = Levenshtein.editops(window, whitelist_seq)
                    differing_base_indices = [op[1] for op in ops if op[0] == "replace"]
                    insertions_or_deletions = [op for op in ops if op[0] in ["insert", "delete"]]

                    if len(insertions_or_deletions) == len(ops):
                        corrected_sequence = "".join(
                            whitelist_seq[op[1]] if op[0] == "replace" else whitelist_seq[op[2]]
                            for op in ops
                        )
                        best_match = {
                            "original_sequence": window,
                            "corrected_sequence": corrected_sequence,
                            "matched_spatial_barcode": whitelist_seq,
                            "match_count": 1,
                            "matched_qualities": quality_windows[i],
                            "levenshtein_distance": lev_dist,
                            "operation_type": "insertion/deletion",
                        }
                        break

                    if differing_base_indices:
                        if all((ord(quality_windows[i][differing_base]) - 33) < 30 for differing_base in differing_base_indices):
                            corrected_sequence = list(window)
                            for differing_base in differing_base_indices:
                                corrected_sequence[differing_base] = whitelist_seq[differing_base]
                            best_match = {
                                "original_sequence": window,
                                "corrected_sequence": "".join(corrected_sequence),
                                "matched_spatial_barcode": whitelist_seq,
                                "match_count": 1,
                                "matched_qualities": quality_windows[i],
                                "levenshtein_distance": lev_dist,
                                "operation_type": "mismatch",
                            }
                            break
                if best_match:
                    break
            if best_match:
                break

        if best_match:
            umi = extract_umi(target_sequence, best_match["original_sequence"])
            updated_row = dict(row)
            updated_row.update(best_match)
            updated_row["UMI"] = umi
            return updated_row
    except Exception as exc:
        print(f"Error while processing row: {exc}")
    return None


def iter_processed_rows(chunk, whitelist):
    # Process one TSV chunk. Keeping the number of submitted futures bounded
    # prevents memory blow-up from queuing all alignment rows at once.
    records = chunk.replace("", float("nan")).to_dict("records")
    if max_workers == 1:
        for row in records:
            result = process_sequence(row, whitelist)
            if result:
                yield result
        return

    pending = set()
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for row in records:
            pending.add(executor.submit(process_sequence, row, whitelist))
            if len(pending) >= max_workers * 4:
                done = {future for future in pending if future.done()}
                for future in done:
                    result = future.result()
                    if result:
                        yield result
                pending -= done
        for future in as_completed(pending):
            result = future.result()
            if result:
                yield result


def append_group(aggregated, columns, row):
    # Aggregate read-level matches by matched spatial barcode. Count columns are
    # summed; other fields are kept as {value}{value} lists.
    barcode = row["matched_spatial_barcode"]
    if barcode not in aggregated:
        aggregated[barcode] = {"match_count": 0, "_values": defaultdict(list)}
    aggregated[barcode]["match_count"] += pd.to_numeric(row.get("match_count", 0), errors="coerce")
    for col in columns:
        if col == "match_count":
            continue
        aggregated[barcode]["_values"][col].append(str(row.get(col, float("nan"))))


def write_aggregated_alignments(aggregated, columns, output_file):
    # Output 1: one row per matched spatial barcode after merging MiXCR fields
    # with spatial metadata.
    rows = []
    for barcode in sorted(aggregated):
        item = aggregated[barcode]
        row = {}
        for col in columns:
            if col == "match_count":
                row[col] = item["match_count"]
            else:
                row[col] = join_braced(item["_values"].get(col, []))
        rows.append(row)
    pd.DataFrame(rows, columns=columns).to_csv(output_file, index=False)


def process_alignment_files(meta_data_file, alignment_files):
    print("[1/8] Matching MiXCR alignments to spatial barcode whitelist and extracting UMI")
    meta_data = pd.read_csv(os.path.join(base_path, meta_data_file))
    whitelist = set(meta_data["spatial_barcode_mixcr"])

    for alignment_file in alignment_files:
        input_file = os.path.join(base_path, f"{alignment_file}-filtered_alignments_poly8A.tsv")
        output_file = os.path.join(base_path, f"{alignment_file}_updated_alignments_umi_poly8A3.csv")
        print(f"Input:  {input_file}")
        print(f"Output: {output_file}")
        aggregated = {}
        output_columns = None
        total_rows = 0
        matched_rows = 0

        for chunk in pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize):
            total_rows += len(chunk)
            updated_records = list(iter_processed_rows(chunk, whitelist))
            matched_rows += len(updated_records)
            if not updated_records:
                print(f"Scanned {total_rows} rows; no matches in current chunk")
                continue

            updated_df = pd.DataFrame(updated_records)
            if "matched_spatial_barcode" not in updated_df.columns:
                print("Error: 'matched_spatial_barcode' column is missing from the updated dataframe.")
                continue

            merged = updated_df.merge(
                meta_data,
                left_on="matched_spatial_barcode",
                right_on="spatial_barcode_mixcr",
                how="left",
            )
            if output_columns is None:
                output_columns = list(merged.columns)
            for row in merged.to_dict("records"):
                append_group(aggregated, output_columns, row)

            print(f"Scanned {total_rows} rows; matched {matched_rows} rows; barcode groups: {len(aggregated)}")

        if output_columns is None:
            pd.DataFrame().to_csv(output_file, index=False)
            print(f"No matches found; saved empty file: {output_file}")
        else:
            write_aggregated_alignments(aggregated, output_columns, output_file)
            print(f"Processing complete; results saved to '{output_file}'")


def extract_unique_barcode(spatial_barcode_ori):
    # spatial_barcode_ori is a braced list; spatial_bc keeps the first barcode
    # as a compact row index for summary tables.
    barcodes = re.findall(r"\{(.*?)\}", spatial_barcode_ori)
    return barcodes[0] if barcodes else None


def filter_by_levenshtein(df, max_levenshtein_distance):
    # Keep only entries whose barcode correction distance is within threshold.
    # The current manuscript workflow uses Lv1.
    filtered_rows = []
    columns = list(df.columns)
    for row in df.to_dict("records"):
        levenshtein_distances = list(map(int, split_braced(row["levenshtein_distance"])))
        filtered_indices = [i for i, distance in enumerate(levenshtein_distances) if distance <= max_levenshtein_distance]
        if not filtered_indices:
            continue

        new_row = {}
        for col in columns:
            value = row[col]
            if col == "levenshtein_distance":
                new_row[col] = join_braced(levenshtein_distances[i] for i in filtered_indices)
            elif isinstance(value, str) and "{" in value and "}" in value:
                entries = split_braced(value)
                new_row[col] = join_braced(entries[i] for i in filtered_indices)
            else:
                new_row[col] = value
        new_row["match_count"] = len(filtered_indices)
        new_row["UMI_count"] = len(set(split_braced(new_row["UMI"])))
        new_row["aaSeqCDR3_count"] = len(set(split_braced(new_row["aaSeqCDR3"])))
        new_row["cloneid_count"] = len(set(split_braced(new_row["cloneId"])))
        filtered_rows.append(new_row)
    return pd.DataFrame(filtered_rows)


def add_levenshtein_outputs(input_files):
    print("[2/8] Applying Levenshtein-distance filter and writing Lv summaries")
    filtered_files = {}
    for input_file in input_files:
        df = pd.read_csv(input_file)
        df["spatial_bc"] = df["spatial_barcode_ori"].apply(extract_unique_barcode)
        df.to_csv(input_file, index=False)
        print(f"Added 'spatial_bc' column and saved updated file: {input_file}")

        for max_dist in max_levenshtein_distances:
            filtered_df = filter_by_levenshtein(df, max_dist)
            output_file = f'{input_file.split(".")[0]}_metadata_with_levenshtein_region_{max_dist}.csv'
            print(f"Output: {output_file}")
            filtered_df.to_csv(output_file, index=False)
            filtered_files[max_dist] = output_file

        reference_file = f'{input_file.split(".")[0]}_metadata_with_levenshtein_region_1.csv'
        reference_df = pd.read_csv(reference_file).set_index("spatial_bc")
        merged_df = reference_df[["UMI_count", "match_count", "aaSeqCDR3_count", "cloneid_count"]].copy()
        merged_df.rename(
            columns={
                "UMI_count": "UMI_count_Lv1",
                "match_count": "match_count_Lv1",
                "aaSeqCDR3_count": "aaSeqCDR3_count_Lv1",
                "cloneid_count": "cloneid_count_Lv1",
            },
            inplace=True,
        )

        for max_dist in max_levenshtein_distances:
            if max_dist == 1:
                continue
            dist_df = pd.read_csv(filtered_files[max_dist]).set_index("spatial_bc")
            merged_df[f"UMI_count_Lv{max_dist}"] = dist_df["UMI_count"]
            merged_df[f"match_count_Lv{max_dist}"] = dist_df["match_count"]
            merged_df[f"aaSeqCDR3_count_Lv{max_dist}"] = dist_df["aaSeqCDR3_count"]
            merged_df[f"cloneid_count_Lv{max_dist}"] = dist_df["cloneid_count"]

        output_file = f'{input_file.split(".")[0]}_TRB_lv_UMI_reads_counts.csv'
        merged_df.to_csv(output_file)
        print(f"Merged data saved to {output_file}")


def calculate_umi_counts(df):
    # Count UMI support within each spatial barcode. These private counts are
    # later used to keep only UMIs observed at least umi_threshold times.
    umi_counts = defaultdict(Counter)
    for spatial_bc, group in df.groupby("spatial_bc"):
        all_umis = re.findall(r"{([^}]+)}", "".join(group["UMI"]))
        umi_counts[spatial_bc] = Counter(all_umis)
    return umi_counts


def add_umi_private_counts(input_files):
    print("[3/8] Counting UMI support within each spatial barcode")
    for umi_threshold in umi_thresholds:
        for max_dist in max_levenshtein_distances:
            for input_file in input_files:
                file_path = f'{input_file.split(".")[0]}_metadata_with_levenshtein_region_{max_dist}.csv'
                output_file_path = f'{input_file.split(".")[0]}_metadata_with_levenshtein_region_umi{umi_threshold}_{max_dist}.csv'
                df = pd.read_csv(file_path)
                umi_counts_by_spatial_bc = calculate_umi_counts(df)

                def replace_counts(row):
                    counts = [umi_counts_by_spatial_bc[row["spatial_bc"]].get(umi, 0) for umi in re.findall(r"{([^}]+)}", row["UMI"])]
                    return join_braced(counts)

                df["UMI_counts_private"] = df.apply(replace_counts, axis=1)
                print(f"Output: {output_file_path}")
                df.to_csv(output_file_path, index=False)


def filter_by_umi_count(df, umi_threshold):
    # Keep read-level entries supported by at least umi_threshold observations
    # of the same UMI within the same spatial barcode.
    filtered_rows = []
    columns = list(df.columns)
    for row in df.to_dict("records"):
        umi_counts_private = list(map(int, split_braced(row["UMI_counts_private"])))
        filtered_indices = [i for i, count in enumerate(umi_counts_private) if count >= umi_threshold]
        if not filtered_indices:
            continue

        new_row = {}
        for col in columns:
            value = row[col]
            if col == "UMI_counts_private":
                new_row[col] = join_braced(umi_counts_private[i] for i in filtered_indices)
            elif isinstance(value, str) and "{" in value and "}" in value:
                entries = split_braced(value)
                new_row[col] = join_braced(entries[i] for i in filtered_indices)
            else:
                new_row[col] = value
        new_row["match_count"] = len(filtered_indices)
        new_row["UMI_count"] = len(set(split_braced(new_row["UMI"])))
        new_row["aaSeqCDR3_count"] = len(set(split_braced(new_row["aaSeqCDR3"])))
        new_row["cloneid_count"] = len(set(split_braced(new_row["cloneId"])))
        filtered_rows.append(new_row)
    return pd.DataFrame(filtered_rows)


def make_filtered_umi_files(input_files):
    print("[4/8] Filtering entries by UMI support threshold")
    for umi_threshold in umi_thresholds:
        for max_dist in max_levenshtein_distances:
            for input_file in input_files:
                input_file_path = f'{input_file.split(".")[0]}_metadata_with_levenshtein_region_umi{umi_threshold}_{max_dist}.csv'
                output_file = f'{input_file.split(".")[0]}_filtered_metadata_with_levenshtein_region_umi{umi_threshold}_{max_dist}.csv'
                df = pd.read_csv(input_file_path)
                print(f"Output: {output_file}")
                filter_by_umi_count(df, umi_threshold).to_csv(output_file, index=False)


def process_filtered_file(input_file_path):
    # Expand braced read-level fields so UMI and cloneId can be resolved. For
    # each spatial barcode + UMI pair, keep the cloneId with the largest support,
    # then collapse back to one row per spatial barcode.
    df = pd.read_csv(input_file_path)
    required_columns = ["UMI_counts_private", "cloneId", "aaSeqCDR3", "spatial_barcode_ori"]
    for col in required_columns:
        if col not in df.columns:
            raise KeyError(f"Error: '{col}' column is missing")

    expanded_columns = {}
    for col in df.columns:
        if df[col].dtype == object and len(df) and "{" in str(df[col].iloc[0]):
            expanded_columns[col] = df[col].str.strip("{}").str.split("}{", expand=True).stack().reset_index(level=1, drop=True)
    expanded_df = pd.concat(expanded_columns, axis=1)
    expanded_df.reset_index(inplace=True)
    expanded_df["UMI_counts_private"] = pd.to_numeric(expanded_df["UMI_counts_private"], errors="coerce")
    expanded_df = expanded_df.dropna(subset=["UMI_counts_private"])
    expanded_df["UMI_counts_private"] = expanded_df["UMI_counts_private"].astype(int)

    most_by_clone = expanded_df.groupby("cloneId")["aaSeqCDR3"].agg(lambda s: Counter(s).most_common(1)[0][0])
    expanded_df["aaSeqCDR3_most"] = expanded_df["cloneId"].map(most_by_clone)

    grouped = expanded_df.groupby(["spatial_barcode_ori", "UMI", "cloneId"]).agg({"UMI_counts_private": "sum"}).reset_index()
    idx = grouped.groupby(["spatial_barcode_ori", "UMI"])["UMI_counts_private"].idxmax()
    max_counts = grouped.loc[idx]
    filtered_df = expanded_df.merge(max_counts[["spatial_barcode_ori", "UMI", "cloneId"]], on=["spatial_barcode_ori", "UMI", "cloneId"])

    rows = []
    for name, group in filtered_df.groupby("spatial_barcode_ori"):
        combined_info = {"spatial_barcode_ori": name}
        for col in group.columns:
            if col != "spatial_barcode_ori":
                combined_info[col] = join_braced(group[col].astype(str))
        rows.append(combined_info)

    final_df = pd.DataFrame(rows)
    final_df["match_count"] = final_df["UMI"].apply(lambda value: len(split_braced(value)))
    final_df["UMI_count"] = final_df["UMI"].apply(lambda value: len(set(split_braced(value))))
    final_df["aaSeqCDR3_count"] = final_df["aaSeqCDR3"].apply(lambda value: len(set(split_braced(value))))
    final_df["cloneid_count"] = final_df["cloneId"].apply(lambda value: len(set(split_braced(value))))
    final_df["aaSeqCDR3_most_count"] = final_df["aaSeqCDR3_most"].apply(lambda value: len(set(split_braced(value))))
    return final_df


def make_cloneid_files(input_files):
    print("[5/8] Resolving cloneId per spatial barcode and UMI")
    for umi_threshold in umi_thresholds:
        for max_dist in max_levenshtein_distances:
            for input_file in input_files:
                input_file_path = f'{input_file.split(".")[0]}_filtered_metadata_with_levenshtein_region_umi{umi_threshold}_{max_dist}.csv'
                output_file = f'{input_file.split(".")[0]}_cloneid_metadata_with_levenshtein_region_umi{umi_threshold}_lv{max_dist}.csv'
                processed_df = process_filtered_file(input_file_path)
                processed_df.to_csv(output_file, index=False)
                print(f"Processed data saved to {output_file}")


def make_filtered_summary_files(input_files):
    print("[6/8] Writing compact barcode-level count summary tables")
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_dist in max_levenshtein_distances:
                file_path = f'{input_file.split(".")[0]}_cloneid_metadata_with_levenshtein_region_umi{umi_threshold}_lv{max_dist}.csv'
                df = pd.read_csv(file_path)
                columns_to_save = [
                    "spatial_barcode_ori",
                    "aaSeqCDR3_most_count",
                    "aaSeqCDR3_count",
                    "cloneid_count",
                    "UMI_count",
                    "match_count",
                ]
                new_file_path = f'filtered_data_{input_file.split(".")[0]}_umi{umi_threshold}_{max_dist}.csv'
                df[columns_to_save].to_csv(new_file_path, index=False, encoding="utf-8")
                print(f"Filtered data saved to {new_file_path}")

    for input_file in input_files:
        merged_df = pd.DataFrame()
        for umi_threshold in umi_thresholds:
            for max_dist in max_levenshtein_distances:
                file_path = f'filtered_data_{input_file.split(".")[0]}_umi{umi_threshold}_{max_dist}.csv'
                df = pd.read_csv(file_path).add_suffix(f"_Lv{max_dist}_UMI{umi_threshold}")
                df.rename(columns={f"spatial_barcode_ori_Lv{max_dist}_UMI{umi_threshold}": "spatial_barcode_ori"}, inplace=True)
                merged_df = df if merged_df.empty else pd.merge(merged_df, df, on="spatial_barcode_ori", how="outer")
        final_output_file = f'merged_filtered_data_{input_file.split(".")[0]}.csv'
        merged_df.to_csv(final_output_file, index=False, encoding="utf-8")
        print(f"Final merged data saved to {final_output_file}")


def generate_cdr3_count_matrix(filtered_df):
    # Matrix A: spatial barcode x aaSeqCDR3_most, counting CDR3 occurrences
    # after Lv/UMI/cloneId filtering.
    matrix_dict = {}
    for row in filtered_df.to_dict("records"):
        if is_missing(row.get("aaSeqCDR3_most")) or is_missing(row.get("spatial_barcode_ori")):
            continue
        spatial_barcode = row["spatial_barcode_ori"]
        bucket = matrix_dict.setdefault(spatial_barcode, {})
        for aaSeqCDR3_most in split_braced(row["aaSeqCDR3_most"]):
            bucket[aaSeqCDR3_most] = bucket.get(aaSeqCDR3_most, 0) + 1
    return pd.DataFrame.from_dict(matrix_dict, orient="index").fillna(0)


def add_spatial_barcode_ori_dub(filtered_df):
    # Duplicate each spatial barcode to align one-to-one with the braced UMI
    # list; this lets the UMI matrix count unique UMIs per barcode/CDR3 pair.
    repeated = []
    for row in filtered_df.to_dict("records"):
        if is_missing(row.get("UMI")):
            repeated.append("")
            continue
        umi_list = split_braced(row["UMI"])
        repeated.append("".join("{" + row["spatial_barcode_ori"] + "}" for _ in umi_list))
    filtered_df["spatial_barcode_ori_dub"] = repeated
    return filtered_df


def generate_umi_count_matrix(filtered_df):
    # Matrix B: spatial barcode x aaSeqCDR3_most, counting unique UMI support.
    matrix_dict = {}
    for row in filtered_df.to_dict("records"):
        if is_missing(row.get("aaSeqCDR3_most")) or is_missing(row.get("spatial_barcode_ori_dub")) or is_missing(row.get("UMI")):
            continue
        for spatial_barcode, aaSeqCDR3_most, umi in zip(
            split_braced(row["spatial_barcode_ori_dub"]),
            split_braced(row["aaSeqCDR3_most"]),
            split_braced(row["UMI"]),
        ):
            matrix_dict.setdefault(spatial_barcode, {}).setdefault(aaSeqCDR3_most, set()).add(umi)
    counts = {key: {sub_key: len(sub_val) for sub_key, sub_val in val.items()} for key, val in matrix_dict.items()}
    return pd.DataFrame.from_dict(counts, orient="index").fillna(0)


def write_matrix_summaries(input_files):
    print("[7/8] Generating aaSeqCDR3 count and UMI count matrices")
    summary_data = []
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_dist in max_levenshtein_distances:
                file_path = f'{input_file.split(".")[0]}_cloneid_metadata_with_levenshtein_region_umi{umi_threshold}_lv{max_dist}.csv'
                filtered_df = pd.read_csv(file_path)
                matrix_df = generate_cdr3_count_matrix(filtered_df)
                output_matrix_file = f'{input_file.split(".")[0]}_aaSeqCDR3_ori_count_matrix_umi{umi_threshold}_levenshtein_{max_dist}.csv'
                matrix_df.to_csv(output_matrix_file)
                summary_data.append({
                    "file_name": f'{input_file.split(".")[0]}_umi{umi_threshold}_levenshtein_{max_dist}',
                    "aaSeqCDR3_count": matrix_df.shape[1],
                    "spatial_barcode_count": matrix_df.shape[0],
                })
    summary_df = pd.DataFrame(summary_data)
    print(summary_df)
    summary_df.to_csv(f'{input_files[-1].split(".")[0]}_unique_aaSeqCDR3_count.csv', index=False, encoding="utf-8")

    summary_data = []
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_dist in max_levenshtein_distances:
                file_path = f'{input_file.split(".")[0]}_cloneid_metadata_with_levenshtein_region_umi{umi_threshold}_lv{max_dist}.csv'
                filtered_df = add_spatial_barcode_ori_dub(pd.read_csv(file_path))
                matrix_df = generate_umi_count_matrix(filtered_df)
                output_matrix_file = f'{input_file.split(".")[0]}_aaSeqCDR3_UMI_count_matrix_umi{umi_threshold}_levenshtein_{max_dist}.csv'
                matrix_df.to_csv(output_matrix_file)
                summary_data.append({
                    "file_name": f'{input_file.split(".")[0]}_umi{umi_threshold}_levenshtein_{max_dist}',
                    "aaSeqCDR3_count": matrix_df.shape[1],
                    "spatial_barcode_count": matrix_df.shape[0],
                })
    summary_df = pd.DataFrame(summary_data)
    print(summary_df)
    summary_df.to_csv(f'{input_files[-1].split(".")[0]}_unique_UMI_aaSeqCDR3_count.csv', index=False, encoding="utf-8")


def calculate_shannon_diversity(matrix_df):
    # Shannon diversity is calculated across CDR3 columns for each spatial
    # barcode using the UMI-count matrix.
    try:
        from skbio.diversity.alpha import shannon
        return matrix_df.apply(lambda row: shannon(row.values), axis=1)
    except ImportError:
        def shannon_from_counts(row):
            values = [float(value) for value in row.values if float(value) > 0]
            total = sum(values)
            return -sum((value / total) * math.log(value / total) for value in values) if total else 0.0
        return matrix_df.apply(shannon_from_counts, axis=1)


def write_shannon(input_files):
    print("[8/8] Calculating Shannon diversity and VDJC/mutation summaries")
    all_shannon_diversities = pd.DataFrame()
    last_input_file = input_files[-1]
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            last_input_file = input_file
            for max_dist in max_levenshtein_distances:
                file_path = f'{input_file.split(".")[0]}_aaSeqCDR3_UMI_count_matrix_umi{umi_threshold}_levenshtein_{max_dist}.csv'
                matrix_df = pd.read_csv(file_path, index_col=0)
                shannon_diversities = calculate_shannon_diversity(matrix_df)
                column_name = f"Shannon_count_Lv{max_dist}_UMI{umi_threshold}"
                if all_shannon_diversities.empty:
                    all_shannon_diversities = pd.DataFrame(shannon_diversities, columns=[column_name])
                else:
                    all_shannon_diversities[column_name] = shannon_diversities
    output_file_path = f'{last_input_file.split(".")[0]}_merged_shannon_diversity.csv'
    all_shannon_diversities.to_csv(output_file_path)
    print("Shannon diversity index calculation and merging completed.")


columns_to_process = [
    "bestVGene", "bestDGene", "bestJGene", "bestCGene", "nSeqCDR3",
    "allVAlignments", "allDAlignments", "allJAlignments", "allCAlignments",
    "UMI",
]


def parse_curly_braces(column_value):
    # Parse braced lists generated earlier, e.g. {TRBV1}{TRBV2}.
    if isinstance(column_value, str):
        return split_braced(column_value)
    return []


def update_genes(group):
    # For each spatial barcode, assign one representative V/D/J/C/nSeqCDR3
    # value per UMI by majority vote among reads sharing that UMI.
    umi_values = group["UMI_parsed"].iloc[0]
    updated_genes = {col: [""] * len(umi_values) for col in columns_to_process if col != "UMI"}

    umi_to_indices = defaultdict(list)
    for idx, umi in enumerate(umi_values):
        umi_to_indices[umi].append(idx)

    for indices in umi_to_indices.values():
        for col in updated_genes:
            gene_values = group[f"{col}_parsed"].iloc[0]
            associated_values = [gene_values[i] for i in indices if i < len(gene_values)]
            most_common_gene = Counter(associated_values).most_common(1)[0][0]
            for i in indices:
                updated_genes[col][i] = most_common_gene

    for col, values in updated_genes.items():
        group[col] = join_braced(values)
    return group


def calculate_unique_umi_count(df, gene_column, umi_column, category):
    # Gene usage frequency is counted by unique UMI, not by raw read count.
    gene_umi_dict = {}
    for row in df.to_dict("records"):
        if is_missing(row.get(gene_column)) or is_missing(row.get(umi_column)):
            continue
        for gene, umi in zip(split_braced(row[gene_column]), split_braced(row[umi_column])):
            gene_umi_dict.setdefault(gene, set()).add(umi)
    frequency_df = pd.DataFrame([(gene, len(umis)) for gene, umis in gene_umi_dict.items()], columns=["Gene", "Frequency"])
    frequency_df["Category"] = category
    return frequency_df


def calculate_vj_paired_unique_umi_count(df):
    # V-J pair frequency is also counted by unique UMI support.
    vj_umi_dict = {}
    for row in df.to_dict("records"):
        if is_missing(row.get("bestVGene")) or is_missing(row.get("bestJGene")) or is_missing(row.get("UMI")):
            continue
        v_genes = split_braced(row["bestVGene"])
        j_genes = split_braced(row["bestJGene"])
        umi_list = split_braced(row["UMI"])
        if len(v_genes) == len(j_genes) == len(umi_list):
            for v_gene, j_gene, umi in zip(v_genes, j_genes, umi_list):
                vj_umi_dict.setdefault(f"{v_gene}_{j_gene}", set()).add(umi)
    frequency_df = pd.DataFrame([(gene, len(umis)) for gene, umis in vj_umi_dict.items()], columns=["Gene", "Frequency"])
    frequency_df["Category"] = "V-J Paired Gene"
    return frequency_df


def calculate_mutation_rate_per_segment(segment):
    # MiXCR alignment strings encode segment coordinates and mutation details.
    # Mutation rate = mutation_count / aligned_segment_length.
    try:
        parts = segment.split("|")
        start, end = int(parts[0]), int(parts[1])
        sequence_length = end - start + 1
        if len(parts) > 5:
            numbers = re.findall(r"\d+", parts[5])
            mutation_count = len(numbers)
            return f"{mutation_count / sequence_length:.4f}", sequence_length, mutation_count
    except Exception:
        return "nan", 0, 0


def calculate_mutation_rate_and_details(value):
    # Return parallel braced lists: mutation rate, aligned length, and mutation
    # count for each V/D/J/C alignment entry.
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
    return join_braced(mutation_rates), join_braced(sequence_lengths), join_braced(mutation_counts)


def expand_rows(df):
    # Expand one barcode-level row back into read/UMI-level rows so regional
    # mutation rates can be calculated per clonotype/UMI observation.
    expanded_rows = []
    for row in df.to_dict("records"):
        spatial_barcode = row["spatial_barcode_ori"]
        lists = {
            "aaSeqCDR3_most": split_braced(row["aaSeqCDR3_most"]),
            "nSeqCDR3": split_braced(row["nSeqCDR3"]),
            "bestVGene": split_braced(row["bestVGene"]),
            "bestDGene": split_braced(row["bestDGene"]),
            "bestJGene": split_braced(row["bestJGene"]),
            "bestCGene": split_braced(row["bestCGene"]),
            "UMI": split_braced(row["UMI"]),
            "allVAlignments_mutation_rate": split_braced(row["allVAlignments_mutation_rate"]),
            "allDAlignments_mutation_rate": split_braced(row["allDAlignments_mutation_rate"]),
            "allJAlignments_mutation_rate": split_braced(row["allJAlignments_mutation_rate"]),
            "allCAlignments_mutation_rate": split_braced(row["allCAlignments_mutation_rate"]),
            "allVAlignments_sequence_length": split_braced(row["allVAlignments_sequence_length"]),
            "allDAlignments_sequence_length": split_braced(row["allDAlignments_sequence_length"]),
            "allJAlignments_sequence_length": split_braced(row["allJAlignments_sequence_length"]),
            "allCAlignments_sequence_length": split_braced(row["allCAlignments_sequence_length"]),
            "allVAlignments_mutation_count": split_braced(row["allVAlignments_mutation_count"]),
            "allDAlignments_mutation_count": split_braced(row["allDAlignments_mutation_count"]),
            "allJAlignments_mutation_count": split_braced(row["allJAlignments_mutation_count"]),
            "allCAlignments_mutation_count": split_braced(row["allCAlignments_mutation_count"]),
        }
        for i in range(len(lists["aaSeqCDR3_most"])):
            expanded_rows.append({"spatial_barcode_ori": spatial_barcode, **{col: values[i] for col, values in lists.items()}})
    return pd.DataFrame(expanded_rows)


def calculate_region_mutation_rates(expanded_df):
    # Combine segment-level mutation counts into VDJ, VDJC, and VJ rates.
    def mutation_rate(row, regions):
        mutation_count = sum(int(row[f"all{region}Alignments_mutation_count"]) for region in regions)
        total_length = sum(int(row[f"all{region}Alignments_sequence_length"]) for region in regions)
        return mutation_count / total_length if total_length > 0 else 0

    expanded_df["VDJ_mutation_rate"] = expanded_df.apply(lambda row: mutation_rate(row, ["V", "D", "J"]), axis=1)
    expanded_df["VDJC_mutation_rate"] = expanded_df.apply(lambda row: mutation_rate(row, ["V", "D", "J", "C"]), axis=1)
    expanded_df["VJ_mutation_rate"] = expanded_df.apply(lambda row: mutation_rate(row, ["V", "J"]), axis=1)


def write_vdjc_outputs(input_files):
    # Outputs:
    #   *_VDJC_umi*_lv*.csv
    #   *_VDJC_VJ_frequency_umi*_lv*.csv
    #   *_with_mutation_rates_umi*_lv*.csv
    #   *_expanded_with_region_rates_umi*_lv*.csv
    for umi_threshold in umi_thresholds:
        for input_file in input_files:
            for max_dist in max_levenshtein_distances:
                file_path = f'{input_file.split(".")[0]}_cloneid_metadata_with_levenshtein_region_umi{umi_threshold}_lv{max_dist}.csv'
                df = pd.read_csv(file_path)

                for col in columns_to_process:
                    df[f"{col}_parsed"] = df[col].apply(parse_curly_braces)

                df = df.groupby("spatial_barcode_ori").apply(update_genes)
                df.drop(columns=[f"{col}_parsed" for col in columns_to_process], inplace=True)

                columns_to_save = [
                    "spatial_barcode_ori", "aaSeqCDR3_most", "aaSeqCDR3_most_count", "aaSeqCDR3_count",
                    "bestVGene", "bestDGene", "bestJGene", "bestCGene",
                    "allVAlignments", "allDAlignments", "allJAlignments", "allCAlignments",
                    "UMI", "UMI_count", "nSeqCDR3",
                ]
                filtered_output_file = f'{input_file.split(".")[0]}_VDJC_umi{umi_threshold}_lv{max_dist}.csv'
                df[columns_to_save].to_csv(filtered_output_file, index=False)

                combined_frequency_df = pd.concat(
                    [
                        calculate_unique_umi_count(df, "bestVGene", "UMI", "VGene"),
                        calculate_unique_umi_count(df, "bestDGene", "UMI", "DGene"),
                        calculate_unique_umi_count(df, "bestJGene", "UMI", "JGene"),
                        calculate_unique_umi_count(df, "bestCGene", "UMI", "CGene"),
                        calculate_vj_paired_unique_umi_count(df),
                    ],
                    ignore_index=True,
                )
                output_file = f'{input_file.split(".")[0]}_VDJC_VJ_frequency_umi{umi_threshold}_lv{max_dist}.csv'
                combined_frequency_df.to_csv(output_file, index=False)
                print(f"VDJC and V-J paired gene frequency calculation completed for {file_path}")

                for col in ["allVAlignments", "allDAlignments", "allJAlignments", "allCAlignments"]:
                    mutation_rates, sequence_lengths, mutation_counts = zip(*df[col].apply(calculate_mutation_rate_and_details))
                    df[f"{col}_mutation_rate"] = mutation_rates
                    df[f"{col}_sequence_length"] = sequence_lengths
                    df[f"{col}_mutation_count"] = mutation_counts

                mutation_rates_output_file = f'{input_file.split(".")[0]}_with_mutation_rates_umi{umi_threshold}_lv{max_dist}.csv'
                df.to_csv(mutation_rates_output_file, index=False)

                expanded_df = expand_rows(df)
                calculate_region_mutation_rates(expanded_df)
                expanded_file_path = f'{input_file.split(".")[0]}_expanded_with_region_rates_umi{umi_threshold}_lv{max_dist}.csv'
                expanded_df.to_csv(expanded_file_path, index=False)
                print(f"Processed {file_path}. Mutation rates saved: {mutation_rates_output_file}, Expanded saved: {expanded_file_path}")


# Execute the workflow in the same order as the original analysis notebook:
# 1) spatial barcode/UMI extraction
# 2) Levenshtein filtering
# 3) UMI support filtering
# 4) cloneId resolution
# 5) summary tables and matrices
# 6) Shannon diversity, VDJC frequency, and mutation-rate exports
for meta_data_file, alignment_files in meta_data_files.items():
    process_alignment_files(meta_data_file, alignment_files)

input_files = [f"{sample_name}_updated_alignments_umi_poly8A3.csv"]
add_levenshtein_outputs(input_files)
add_umi_private_counts(input_files)
make_filtered_umi_files(input_files)
make_cloneid_files(input_files)
make_filtered_summary_files(input_files)
write_matrix_summaries(input_files)
write_shannon(input_files)
write_vdjc_outputs(input_files)
PY

done
