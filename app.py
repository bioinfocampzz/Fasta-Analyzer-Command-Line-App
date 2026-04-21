import argparse
import os
from typing import Dict, List


# =========================
# File Handling Functions
# =========================

def validate_file(file_path: str) -> None:
    """
    Validate file existence and basic properties.
    Raises exceptions if invalid.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    if os.path.getsize(file_path) == 0:
        raise ValueError(f"File is empty: {file_path}")


def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Parse FASTA file and return dictionary {id: sequence}.
    Ensures valid FASTA format.
    """
    validate_file(file_path)

    sequences = {}
    seq_id = None
    seq_chunks = []

    with open(file_path, "r") as file:
        first_line_checked = False

        for line in file:
            line = line.strip()

            if not line:
                continue

            # Validate first meaningful line
            if not first_line_checked:
                if not line.startswith(">"):
                    raise ValueError(f"Invalid FASTA format: {file_path}")
                first_line_checked = True

            if line.startswith(">"):
                # Save previous sequence
                if seq_id:
                    sequences[seq_id] = ''.join(seq_chunks)

                seq_id = line[1:].split()[0]
                if not seq_id:
                    raise ValueError(f"Missing sequence ID in {file_path}")

                seq_chunks = []
            else:
                if seq_id is None:
                    raise ValueError(f"Sequence before header in {file_path}")
                seq_chunks.append(line)

        # Save last sequence
        if seq_id:
            sequences[seq_id] = ''.join(seq_chunks)

    if not sequences:
        raise ValueError(f"No valid sequences found in {file_path}")

    return sequences


# =========================
# Sequence Processing
# =========================

DNA_BASES = {'A', 'T', 'G', 'C'}
PROTEIN_BASES = {
    'A','R','N','D','C','Q','E','G','H','I',
    'L','K','M','F','P','S','T','W','Y','V'
}


def detect_sequence_type(sequence: str) -> str:
    """
    Detect whether sequence is DNA or Protein.
    """
    seq = sequence.upper()

    if set(seq).issubset(DNA_BASES):
        return "DNA"
    elif set(seq).issubset(PROTEIN_BASES):
        return "Protein"
    else:
        return "Invalid"


def validate_sequence(sequence: str, seq_type: str) -> None:
    """
    Validate sequence characters based on type.
    """
    if not sequence:
        raise ValueError("Sequence is empty")

    seq = sequence.upper()

    valid_set = DNA_BASES if seq_type == "DNA" else PROTEIN_BASES

    invalid_chars = set(seq) - valid_set
    if invalid_chars:
        raise ValueError(f"Invalid characters: {', '.join(invalid_chars)}")


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate GC percentage for DNA sequence.
    """
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    return round((gc_count / len(seq)) * 100, 2)


# =========================
# Analysis Logic
# =========================

def analyze_sequence(seq_id: str, sequence: str) -> dict:
    """
    Analyze a single sequence and return result dictionary.
    """
    try:
        seq_type = detect_sequence_type(sequence)

        if seq_type == "Invalid":
            return {
                "id": seq_id,
                "type": "Invalid",
                "length": len(sequence),
                "gc": "N/A"
            }

        validate_sequence(sequence, seq_type)

        gc_value = calculate_gc_content(sequence) if seq_type == "DNA" else "N/A"

        return {
            "id": seq_id,
            "type": seq_type,
            "length": len(sequence),
            "gc": gc_value
        }

    except ValueError as e:
        return {
            "id": seq_id,
            "type": "Error",
            "length": "-",
            "gc": str(e)
        }


def analyze_file(file_path: str) -> List[dict]:
    """
    Analyze all sequences in a FASTA file.
    """
    try:
        sequences = read_fasta(file_path)
        return [analyze_sequence(seq_id, seq) for seq_id, seq in sequences.items()]

    except Exception as e:
        print(f"❌ Error in file '{file_path}': {e}")
        return []


# =========================
# Output Formatting
# =========================

def print_table(results: List[dict], file_name: str) -> None:
    """
    Print formatted table of results.
    """
    if not results:
        print(f"⚠️ No valid results for {file_name}")
        return

    print(f"\n📄 File: {file_name}")
    print("-" * 70)
    print(f"{'Sequence ID':<20}{'Type':<10}{'Length':<10}{'GC %':<10}")
    print("-" * 70)

    for row in results:
        print(f"{row['id']:<20}{row['type']:<10}{row['length']:<10}{row['gc']:<10}")

    print("-" * 70)


# =========================
# CLI Entry Point
# =========================

def parse_arguments():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="FASTA Analyzer Tool (Refactored & Modular)"
    )

    parser.add_argument(
        "fasta_files",
        nargs="+",
        help="One or more FASTA files"
    )

    return parser.parse_args()


def main():
    """
    Main execution function.
    """
    args = parse_arguments()

    for file_path in args.fasta_files:
        print(f"\nProcessing: {file_path}")

        results = analyze_file(file_path)

        if results:
            print_table(results, file_path)
        else:
            print(f"⚠️ Skipping file: {file_path}")


if __name__ == "__main__":
    main()