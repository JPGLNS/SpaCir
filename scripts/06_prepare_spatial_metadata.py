#!/usr/bin/env python3
"""Build SpaCir spatial-barcode metadata from Space Ranger barcodes."""

import argparse
import csv
import gzip
import os
import re
import sys
import tempfile
from pathlib import Path


BARCODE_PATTERN = re.compile(r"^([ACGT]{16})(-\d+)$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def read_barcodes(path):
    records = []
    original_seen = set()
    mixcr_seen = set()

    try:
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
        with handle:
            for line_number, raw_line in enumerate(handle, start=1):
                barcode = raw_line.rstrip("\r\n")
                match = BARCODE_PATTERN.fullmatch(barcode)
                if match is None:
                    raise ValueError(
                        f"invalid barcode at line {line_number}: {barcode!r}; "
                        "expected 16 A/C/G/T bases followed by a numeric suffix"
                    )

                if barcode in original_seen:
                    raise ValueError(
                        f"duplicate original barcode at line {line_number}: {barcode}"
                    )

                mixcr_barcode = reverse_complement(match.group(1))
                if mixcr_barcode in mixcr_seen:
                    raise ValueError(
                        f"duplicate MiXCR barcode at line {line_number}: "
                        f"{mixcr_barcode}"
                    )

                original_seen.add(barcode)
                mixcr_seen.add(mixcr_barcode)
                records.append((barcode, mixcr_barcode))
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read gzip barcode file {path}: {exc}") from exc

    if not records:
        raise ValueError("barcode file contains no records")
    return records


def validate_sample_id(sample_id):
    if not sample_id or sample_id in {".", ".."} or "/" in sample_id or "\\" in sample_id:
        raise ValueError(
            f"invalid sample ID {sample_id!r}: path separators are not allowed"
        )


def write_metadata(records, output_path, force=False):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and not force:
        raise FileExistsError(
            f"output already exists: {output_path}; use --force to replace it"
        )

    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="",
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["spatial_barcode_ori", "spatial_barcode_mixcr"])
            writer.writerows(records)

        if force:
            os.replace(temporary_path, output_path)
        else:
            os.link(temporary_path, output_path)
            temporary_path.unlink()
        temporary_path = None
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Convert Space Ranger barcodes.tsv.gz into SpaCir spatial "
            "barcode metadata."
        )
    )
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--force",
        action="store_true",
        help="replace an existing metadata CSV",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    try:
        validate_sample_id(args.sample_id)
        records = read_barcodes(args.barcodes)
        output_path = (
            args.output_dir / f"{args.sample_id}_meta_data_with_mixcr.csv"
        )
        write_metadata(records, output_path, force=args.force)
    except (FileExistsError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    print(f"Saved spatial metadata: {output_path}")
    print(f"barcodes={len(records)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
