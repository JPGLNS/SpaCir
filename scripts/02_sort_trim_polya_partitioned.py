#!/usr/bin/env python3
"""Sort primer orientations, filter read pairs by polyA signal, and trim read1."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import tempfile
import zlib
from collections import Counter
from contextlib import ExitStack
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq


def load_env(path: Path) -> dict[str, str]:
    env: dict[str, str] = {}
    if not path.exists():
        return env
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        env[key] = value.strip().strip('"')
    return env


def read_samples(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def core_id(seq_id: str) -> str:
    return seq_id.split("/")[0]


def primer_type(record) -> str:
    match = re.search(r"\|PRIMER=([^_\s|]+)", record.description)
    return match.group(1) if match else "unknown"


def max_consecutive_a(sequence: str) -> int:
    runs = re.findall(r"A+", sequence)
    return max(map(len, runs), default=0)


def has_polya_signal(sequence: str, min_run: int, search_bp: int) -> bool:
    return re.search(f"A{{{min_run},}}", sequence[:search_bp]) is not None


def trim_after_last_polya(sequence: str, qualities: list[int], min_run: int, search_bp: int) -> tuple[str, list[int]]:
    matches = list(re.finditer(f"A{{{min_run},}}", sequence[:search_bp]))
    if not matches:
        return sequence, qualities
    trim_pos = matches[-1].end()
    return sequence[trim_pos:], qualities[trim_pos:]


def write_histogram(counts: Counter[int], output_pdf: Path, title: str) -> None:
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(8, 5))
    values = sorted(counts)
    weights = [counts[value] for value in values]
    bins = range(0, 41)
    plt.hist(
        values,
        bins=bins,
        weights=weights,
        color="#4C78A8",
        edgecolor="black",
        align="left",
    )
    plt.title(title)
    plt.xlabel("Maximum consecutive A bases")
    plt.ylabel("Read count")
    plt.xticks(range(0, 41, 5))
    plt.xlim(0, 40)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_pdf)
    plt.close()


def output_paths(output_prefix: Path) -> dict[str, Path]:
    return {
        "rf_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read1.fastq"),
        "rf_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read2.fastq"),
        "fr_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read1.fastq"),
        "fr_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read2.fastq"),
        "both_f": output_prefix.with_name(output_prefix.name + "_both_F.fastq"),
        "both_r": output_prefix.with_name(output_prefix.name + "_both_R.fastq"),
        "no_match": output_prefix.with_name(output_prefix.name + "_no_match.fastq"),
        "rf_poly_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read1_poly8A.fastq"),
        "rf_poly_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read2_poly8A.fastq"),
        "rf_trimmed_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_R_F_read1_poly8A_trimmed.fastq"),
        "fr_poly_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read1_poly8A.fastq"),
        "fr_poly_r2": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read2_poly8A.fastq"),
        "fr_trimmed_r1": output_prefix.with_name(output_prefix.name + "_PRIMER_F_R_read1_poly8A_trimmed.fastq"),
    }

def count_fastq_records(path: Path) -> int:
    line_count = 0
    with path.open("rb") as handle:
        for line_count, _ in enumerate(handle, start=1):
            pass
    if line_count % 4:
        raise ValueError(f"FASTQ is not a four-line record file: {path}")
    return line_count // 4


def fastq_partition(header: bytes, partition_count: int) -> int:
    read_id = header[1:].split(None, 1)[0].split(b"/", 1)[0]
    return zlib.crc32(read_id) % partition_count


def partition_fastq(path: Path, output_dir: Path, partition_count: int) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = [output_dir / f"part{index:04d}.fastq" for index in range(partition_count)]
    with ExitStack() as stack:
        outputs = [stack.enter_context(item.open("wb")) for item in paths]
        with path.open("rb") as source:
            record_count = 0
            while header := source.readline():
                sequence = source.readline()
                separator = source.readline()
                quality = source.readline()
                if not header.startswith(b"@") or not separator.startswith(b"+") or not quality:
                    raise ValueError(f"Invalid or truncated FASTQ record in {path}")
                partition = fastq_partition(header, partition_count)
                outputs[partition].writelines((header, sequence, separator, quality))
                record_count += 1
                if record_count % 1_000_000 == 0:
                    print(f"Partitioned {record_count:,} records from {path.name}", flush=True)
    return paths


def write_polya_pair(record1, record2, label: str, handles, stats, min_run: int, search_bp: int) -> None:
    sequence = str(record1.seq)
    count = max_consecutive_a(sequence)
    stats[label]["all"][count] += 1
    if not has_polya_signal(sequence, min_run, search_bp):
        return

    stats[label]["kept"][count] += 1
    SeqIO.write(record1, handles[f"{label}_poly_r1"], "fastq")
    SeqIO.write(record2, handles[f"{label}_poly_r2"], "fastq")
    qualities = record1.letter_annotations["phred_quality"]
    trimmed_sequence, trimmed_quality = trim_after_last_polya(sequence, qualities, min_run, search_bp)
    record1.letter_annotations = {}
    record1.seq = Seq(trimmed_sequence)
    record1.letter_annotations["phred_quality"] = trimmed_quality
    SeqIO.write(record1, handles[f"{label}_trimmed_r1"], "fastq")


def process_partition(read1_path: Path, read2_path: Path, handles, stats, min_run: int, search_bp: int) -> None:
    # Only one R2 partition is indexed. R1 is streamed and matched records are
    # removed immediately, keeping the peak proportional to one partition.
    with read2_path.open() as read2_handle:
        reads2 = {core_id(record.id): record for record in SeqIO.parse(read2_handle, "fastq")}
    with read1_path.open() as read1_handle:
        for record1 in SeqIO.parse(read1_handle, "fastq"):
            record2 = reads2.pop(core_id(record1.id), None)
            if record2 is None:
                SeqIO.write(record1, handles["no_match"], "fastq")
                continue

            p1 = primer_type(record1)
            p2 = primer_type(record2)
            if p1 == "Fprimer" and p2 == "Rprimer":
                normalized1, normalized2, label = record1, record2, "fr"
            elif p1 == "Rprimer" and p2 == "Fprimer":
                normalized1, normalized2, label = record2, record1, "rf"
            elif p1 == "Fprimer" and p2 == "Fprimer":
                SeqIO.write((record1, record2), handles["both_f"], "fastq")
                continue
            elif p1 == "Rprimer" and p2 == "Rprimer":
                SeqIO.write((record1, record2), handles["both_r"], "fastq")
                continue
            else:
                SeqIO.write((record1, record2), handles["no_match"], "fastq")
                continue

            SeqIO.write(normalized1, handles[f"{label}_r1"], "fastq")
            SeqIO.write(normalized2, handles[f"{label}_r2"], "fastq")
            write_polya_pair(normalized1, normalized2, label, handles, stats, min_run, search_bp)

    for record2 in reads2.values():
        SeqIO.write(record2, handles["no_match"], "fastq")


def process_partitioned_sample(
    read1_path: Path,
    read2_path: Path,
    output_prefix: Path,
    min_run: int,
    search_bp: int,
    target_records: int,
) -> None:
    if target_records <= 0:
        raise ValueError("STAGE2_PARTITION_READS must be greater than zero")
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    read1_count = count_fastq_records(read1_path)
    read2_count = count_fastq_records(read2_path)
    partition_count = max(1, math.ceil(max(read1_count, read2_count) / target_records))
    print(
        f"Stage 2 input: R1={read1_count:,}, R2={read2_count:,}; "
        f"using {partition_count} ID partitions (target {target_records:,} records each)",
        flush=True,
    )

    temp_root = Path(tempfile.mkdtemp(prefix=f".{output_prefix.name}_stage2_", dir=output_prefix.parent))
    try:
        read1_parts = partition_fastq(read1_path, temp_root / "read1", partition_count)
        read2_parts = partition_fastq(read2_path, temp_root / "read2", partition_count)
        paths = output_paths(output_prefix)
        stats = {
            "fr": {"all": Counter(), "kept": Counter()},
            "rf": {"all": Counter(), "kept": Counter()},
        }
        with ExitStack() as stack:
            handles = {name: stack.enter_context(path.open("w")) for name, path in paths.items()}
            for index, (read1_part, read2_part) in enumerate(zip(read1_parts, read2_parts), start=1):
                print(f"Processing ID partition {index}/{partition_count}", flush=True)
                process_partition(read1_part, read2_part, handles, stats, min_run, search_bp)
                read1_part.unlink()
                read2_part.unlink()

        for label in ("fr", "rf"):
            plot_prefix = output_prefix.with_name(f"{output_prefix.name}_{label}")
            read1_name = paths[f"{label}_r1"].name
            write_histogram(
                stats[label]["all"],
                plot_prefix.with_name(plot_prefix.name + "_polya_all.pdf"),
                f"polyA distribution: {read1_name}",
            )
            write_histogram(
                stats[label]["kept"],
                plot_prefix.with_name(plot_prefix.name + "_polya_kept.pdf"),
                f"polyA+ reads: {read1_name}",
            )
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--env", default="config/pipeline.env")
    args = parser.parse_args()

    env = load_env(Path(args.env))
    results_dir = Path(env.get("RESULTS_DIR", "results"))
    preprocess_dir = env.get("PREPROCESS_DIR", "preprocessing")
    min_run = int(env.get("POLYA_MIN_RUN", "8"))
    search_bp = int(env.get("POLYA_SEARCH_BP", "80"))
    target_records = int(env.get("STAGE2_PARTITION_READS", "2000000"))

    for sample in read_samples(Path(args.samples)):
        sample_id = sample["sample_id"]
        sample_dir = results_dir / sample_id / preprocess_dir
        output_prefix = sample_dir / f"{sample_id}_output"
        read1_pass = sample_dir / (Path(sample["read1"]).name.replace(".gz", "").replace(".fastq", "").replace(".fq", "") + "_primers-pass.fastq")
        read2_pass = sample_dir / (Path(sample["read2"]).name.replace(".gz", "").replace(".fastq", "").replace(".fq", "") + "_primers-pass.fastq")

        process_partitioned_sample(
            read1_pass,
            read2_pass,
            output_prefix,
            min_run,
            search_bp,
            target_records,
        )
        print(f"Finished primer sorting and polyA trimming for {sample_id}")


if __name__ == "__main__":
    main()
