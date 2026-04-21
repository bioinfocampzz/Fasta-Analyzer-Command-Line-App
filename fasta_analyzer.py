#!/usr/bin/env python3
"""
fasta_analyzer.py — Production-quality CLI tool for FASTA file analysis.

Supports DNA and protein sequences. Provides per-sequence stats (length,
GC content, amino acid composition), formatted table output, and optional
CSV export.

Usage:
    python fasta_analyzer.py file.fasta
    python fasta_analyzer.py *.fasta --type DNA --csv results.csv
    python fasta_analyzer.py /data/genomes/ --recursive
"""

import argparse
import csv
import os
import sys
import warnings
from collections import Counter
from pathlib import Path
from typing import Generator, Optional

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna", ".faa"}

DNA_VALID_CHARS = set("ATGCN")
PROTEIN_VALID_CHARS = set("ACDEFGHIKLMNPQRSTVWYX")

# Standard one-letter amino acid codes
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

class FastaRecord:
    """Lightweight container for a single FASTA sequence record."""

    __slots__ = ("seq_id", "description", "sequence", "source_file")

    def __init__(
        self,
        seq_id: str,
        description: str,
        sequence: str,
        source_file: str,
    ) -> None:
        self.seq_id = seq_id
        self.description = description
        self.sequence = sequence
        self.source_file = source_file

    def __repr__(self) -> str:
        return f"<FastaRecord id={self.seq_id!r} len={len(self.sequence)}>"


class AnalysisResult:
    """Holds computed statistics for one sequence."""

    __slots__ = (
        "file_name",
        "seq_id",
        "seq_type",
        "length",
        "gc_content",
        "aa_composition",
        "invalid_chars",
    )

    def __init__(
        self,
        file_name: str,
        seq_id: str,
        seq_type: str,
        length: int,
        gc_content: Optional[float],
        aa_composition: Optional[dict],
        invalid_chars: set,
    ) -> None:
        self.file_name = file_name
        self.seq_id = seq_id
        self.seq_type = seq_type
        self.length = length
        self.gc_content = gc_content
        self.aa_composition = aa_composition
        self.invalid_chars = invalid_chars


# ---------------------------------------------------------------------------
# File collection
# ---------------------------------------------------------------------------

def collect_fasta_files(paths: list[str], recursive: bool = False) -> list[Path]:
    """
    Collect all FASTA files from a list of file paths or directories.

    Args:
        paths:     List of file or directory path strings supplied by the user.
        recursive: Whether to scan directories recursively.

    Returns:
        Sorted, deduplicated list of Path objects pointing to FASTA files.
    """
    collected: list[Path] = []

    for raw in paths:
        p = Path(raw)

        if p.is_file():
            if p.suffix.lower() in FASTA_EXTENSIONS:
                collected.append(p)
            else:
                _warn(f"Skipping '{p}': unrecognised extension (expected {FASTA_EXTENSIONS})")

        elif p.is_dir():
            pattern = "**/*" if recursive else "*"
            for candidate in p.glob(pattern):
                if candidate.is_file() and candidate.suffix.lower() in FASTA_EXTENSIONS:
                    collected.append(candidate)

        else:
            _warn(f"Path not found: '{raw}'")

    # Deduplicate while preserving order
    seen: set[Path] = set()
    unique: list[Path] = []
    for f in collected:
        resolved = f.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique.append(f)

    return sorted(unique)


# ---------------------------------------------------------------------------
# FASTA reading (streaming — memory-efficient)
# ---------------------------------------------------------------------------

def read_fasta(file_path: Path) -> Generator[FastaRecord, None, None]:
    """
    Stream-parse a FASTA file, yielding one FastaRecord per sequence.

    - Skips blank lines.
    - Handles multi-line sequences.
    - Raises ValueError on format violations.

    Args:
        file_path: Path to the FASTA file.

    Yields:
        FastaRecord objects.

    Raises:
        ValueError: If the file is empty or a sequence header is missing.
        FileNotFoundError: If the file does not exist.
    """
    current_header: Optional[str] = None
    seq_parts: list[str] = []
    has_any_record = False

    with file_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line_no, raw_line in enumerate(fh, start=1):
            line = raw_line.strip()

            if not line:
                continue  # Skip blank / whitespace-only lines

            if line.startswith(">"):
                # Flush the previous record before starting a new one
                if current_header is not None:
                    yield _build_record(current_header, seq_parts, file_path)
                    seq_parts = []

                current_header = line[1:].strip()  # Strip the leading ">"
                has_any_record = True

            else:
                if current_header is None:
                    raise ValueError(
                        f"[{file_path.name}] Line {line_no}: sequence data found before "
                        "any header line. Not a valid FASTA file."
                    )
                seq_parts.append(line.upper())

    # Flush the last record
    if current_header is not None:
        yield _build_record(current_header, seq_parts, file_path)
    elif not has_any_record:
        raise ValueError(f"[{file_path.name}] File is empty or contains no FASTA records.")


def _build_record(header: str, seq_parts: list[str], file_path: Path) -> FastaRecord:
    """
    Assemble a FastaRecord from a collected header and sequence fragments.

    Args:
        header:    Full header string (without ">").
        seq_parts: List of sequence line strings.
        file_path: Source file (for attribution).

    Returns:
        A completed FastaRecord.

    Raises:
        ValueError: When the sequence body is empty.
    """
    tokens = header.split(maxsplit=1)
    seq_id = tokens[0] if tokens else "(unnamed)"
    description = tokens[1] if len(tokens) > 1 else ""
    sequence = "".join(seq_parts)

    if not sequence:
        raise ValueError(
            f"[{file_path.name}] Header '>{header}' has no associated sequence."
        )

    return FastaRecord(
        seq_id=seq_id,
        description=description,
        sequence=sequence,
        source_file=str(file_path),
    )


# ---------------------------------------------------------------------------
# Sequence utilities
# ---------------------------------------------------------------------------

def clean_sequence(sequence: str) -> str:
    """
    Remove whitespace and convert sequence to uppercase.

    Args:
        sequence: Raw sequence string.

    Returns:
        Cleaned, uppercased sequence with no whitespace.
    """
    return "".join(sequence.split()).upper()


def detect_sequence_type(sequence: str) -> str:
    """
    Heuristically determine whether a sequence is DNA or Protein.

    A sequence is classified as DNA when ≥ 85 % of its characters belong to
    {A, T, G, C, N}.  This threshold accommodates minor contamination while
    avoiding false positives for protein sequences that happen to share single-
    letter codes with nucleotides.

    Args:
        sequence: Cleaned, uppercased sequence string.

    Returns:
        "DNA" or "PROTEIN".
    """
    if not sequence:
        return "UNKNOWN"

    dna_chars = sum(1 for c in sequence if c in DNA_VALID_CHARS)
    dna_fraction = dna_chars / len(sequence)

    return "DNA" if dna_fraction >= 0.85 else "PROTEIN"


def find_invalid_chars(sequence: str, seq_type: str) -> set[str]:
    """
    Identify characters in the sequence that are invalid for the given type.

    Args:
        sequence: Cleaned sequence string.
        seq_type: "DNA" or "PROTEIN".

    Returns:
        Set of invalid character strings (may be empty).
    """
    valid = DNA_VALID_CHARS if seq_type == "DNA" else PROTEIN_VALID_CHARS
    return {c for c in sequence if c not in valid}


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def calculate_gc(sequence: str) -> float:
    """
    Compute GC content as a percentage.

    Args:
        sequence: Cleaned DNA sequence (uppercase).

    Returns:
        GC percentage in range [0.0, 100.0].
        Returns 0.0 for an empty sequence.
    """
    if not sequence:
        return 0.0

    gc_count = sequence.count("G") + sequence.count("C")
    return round(gc_count / len(sequence) * 100, 2)


def amino_acid_composition(sequence: str) -> dict[str, int]:
    """
    Count occurrences of each standard amino acid in a protein sequence.

    Non-standard characters are counted under the key "OTHER".

    Args:
        sequence: Cleaned protein sequence (uppercase).

    Returns:
        Dictionary mapping single-letter amino acid codes to counts.
    """
    counts = Counter(sequence)
    composition: dict[str, int] = {}

    for aa in AMINO_ACIDS:
        composition[aa] = counts.get(aa, 0)

    # Aggregate non-standard residues
    other = sum(v for k, v in counts.items() if k not in set(AMINO_ACIDS))
    if other:
        composition["OTHER"] = other

    return composition


def generate_stats(record: FastaRecord, seq_type: str) -> AnalysisResult:
    """
    Run the full analysis pipeline for one FASTA record.

    Args:
        record:   A parsed FastaRecord.
        seq_type: "DNA" or "PROTEIN" (already resolved or auto-detected).

    Returns:
        An AnalysisResult containing all computed statistics.
    """
    sequence = clean_sequence(record.sequence)
    invalid = find_invalid_chars(sequence, seq_type)

    if invalid:
        _warn(
            f"Sequence '{record.seq_id}' in '{Path(record.source_file).name}' "
            f"contains invalid {seq_type} characters: {sorted(invalid)}"
        )

    gc = calculate_gc(sequence) if seq_type == "DNA" else None
    aa_comp = amino_acid_composition(sequence) if seq_type == "PROTEIN" else None

    return AnalysisResult(
        file_name=Path(record.source_file).name,
        seq_id=record.seq_id,
        seq_type=seq_type,
        length=len(sequence),
        gc_content=gc,
        aa_composition=aa_comp,
        invalid_chars=invalid,
    )


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _warn(message: str) -> None:
    """Print a yellow-prefixed warning to stderr."""
    print(f"\033[33m[WARNING]\033[0m {message}", file=sys.stderr)


def _error(message: str) -> None:
    """Print a red-prefixed error to stderr."""
    print(f"\033[31m[ERROR]\033[0m {message}", file=sys.stderr)


def print_results_table(results: list[AnalysisResult]) -> None:
    """
    Render analysis results as a formatted, aligned ASCII table.

    Args:
        results: List of AnalysisResult objects.
    """
    if not results:
        print("No results to display.")
        return

    # Column headers
    headers = ["File Name", "Sequence ID", "Type", "Length", "GC Content (%)"]

    # Build rows
    rows: list[list[str]] = []
    for r in results:
        gc_str = f"{r.gc_content:.2f}" if r.gc_content is not None else "N/A"
        rows.append([r.file_name, r.seq_id, r.seq_type, str(r.length), gc_str])

    # Compute column widths
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(cell))

    sep = "+-" + "-+-".join("-" * w for w in col_widths) + "-+"
    header_row = "| " + " | ".join(h.ljust(col_widths[i]) for i, h in enumerate(headers)) + " |"

    print()
    print(sep)
    print(header_row)
    print(sep)
    for row in rows:
        print("| " + " | ".join(cell.ljust(col_widths[i]) for i, cell in enumerate(row)) + " |")
    print(sep)
    print()


def print_summary(results: list[AnalysisResult]) -> None:
    """
    Print aggregate summary statistics to stdout.

    Args:
        results: List of AnalysisResult objects.
    """
    total = len(results)
    dna_count = sum(1 for r in results if r.seq_type == "DNA")
    protein_count = sum(1 for r in results if r.seq_type == "PROTEIN")

    print("=" * 50)
    print("  SUMMARY")
    print("=" * 50)
    print(f"  Total sequences analysed : {total}")
    print(f"  DNA sequences            : {dna_count}")
    print(f"  Protein sequences        : {protein_count}")

    if dna_count:
        dna_results = [r for r in results if r.seq_type == "DNA"]
        avg_gc = sum(r.gc_content for r in dna_results) / dna_count
        print(f"  Average GC content (DNA) : {avg_gc:.2f}%")

    print("=" * 50)
    print()


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def export_csv(results: list[AnalysisResult], output_path: str) -> None:
    """
    Write all results to a CSV file.

    Protein amino acid composition is serialised as a compact key:value string
    so that every record fits in a single CSV row.

    Args:
        results:     List of AnalysisResult objects.
        output_path: Destination file path.
    """
    fieldnames = [
        "file_name",
        "seq_id",
        "seq_type",
        "length",
        "gc_content",
        "aa_composition",
        "invalid_chars",
    ]

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()

        for r in results:
            aa_str = (
                ";".join(f"{k}:{v}" for k, v in r.aa_composition.items())
                if r.aa_composition
                else ""
            )
            writer.writerow(
                {
                    "file_name": r.file_name,
                    "seq_id": r.seq_id,
                    "seq_type": r.seq_type,
                    "length": r.length,
                    "gc_content": r.gc_content if r.gc_content is not None else "",
                    "aa_composition": aa_str,
                    "invalid_chars": "".join(sorted(r.invalid_chars)),
                }
            )

    print(f"Results exported to: {output_path}")


# ---------------------------------------------------------------------------
# CLI argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """
    Construct and return the argument parser.

    Returns:
        Configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(
        prog="fasta_analyzer",
        description=(
            "Analyse FASTA files for DNA/Protein sequences.\n"
            "Computes GC content (DNA) or amino acid composition (Protein),\n"
            "prints a formatted table, and optionally exports to CSV."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python fasta_analyzer.py genome.fasta\n"
            "  python fasta_analyzer.py *.fa --type DNA --csv results.csv\n"
            "  python fasta_analyzer.py /data/seqs/ --recursive\n"
        ),
    )

    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="INPUT",
        help="One or more FASTA files or directories.",
    )
    parser.add_argument(
        "--type",
        choices=["DNA", "PROTEIN"],
        default=None,
        metavar="TYPE",
        help="Force sequence type (DNA or PROTEIN). Auto-detected if omitted.",
    )
    parser.add_argument(
        "--csv",
        metavar="FILE",
        default=None,
        help="Export results to this CSV file.",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Scan directories recursively for FASTA files.",
    )
    parser.add_argument(
        "--no-warnings",
        action="store_true",
        help="Suppress warning messages.",
    )

    return parser


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """
    Application entry point.

    Orchestrates argument parsing, file collection, sequence analysis,
    table rendering, and optional CSV export.
    """
    parser = build_parser()
    args = parser.parse_args()

    # Optionally suppress warnings by monkey-patching _warn
    if args.no_warnings:
        global _warn
        _warn = lambda _: None  # noqa: E731

    # ---- Collect files -------------------------------------------------------
    fasta_files = collect_fasta_files(args.inputs, recursive=args.recursive)

    if not fasta_files:
        _error("No valid FASTA files found. Check your input paths and extensions.")
        sys.exit(1)

    print(f"\nFound {len(fasta_files)} FASTA file(s) to process.\n")

    # ---- Process sequences ---------------------------------------------------
    all_results: list[AnalysisResult] = []
    files_with_errors: list[str] = []

    for fasta_path in fasta_files:
        print(f"Processing: {fasta_path}")
        try:
            for record in read_fasta(fasta_path):
                # Determine sequence type
                seq_type = args.type if args.type else detect_sequence_type(record.sequence)

                result = generate_stats(record, seq_type)
                all_results.append(result)

        except (ValueError, UnicodeDecodeError) as exc:
            _error(str(exc))
            files_with_errors.append(str(fasta_path))

        except FileNotFoundError:
            _error(f"File disappeared during processing: '{fasta_path}'")
            files_with_errors.append(str(fasta_path))

        except Exception as exc:  # pylint: disable=broad-except
            _error(f"Unexpected error processing '{fasta_path}': {exc}")
            files_with_errors.append(str(fasta_path))

    # ---- Guard: nothing successfully analysed --------------------------------
    if not all_results:
        _error("No sequences could be analysed. Exiting.")
        sys.exit(1)

    # ---- Render output -------------------------------------------------------
    print_results_table(all_results)
    print_summary(all_results)

    if files_with_errors:
        print(f"  Files skipped due to errors: {len(files_with_errors)}")
        for f in files_with_errors:
            print(f"    - {f}")
        print()

    # ---- CSV export ----------------------------------------------------------
    if args.csv:
        try:
            export_csv(all_results, args.csv)
        except OSError as exc:
            _error(f"Could not write CSV file '{args.csv}': {exc}")
            sys.exit(1)


if __name__ == "__main__":
    main()