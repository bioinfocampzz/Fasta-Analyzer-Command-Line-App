# FASTA Analyzer

A production-quality command-line tool for analyzing FASTA files containing DNA, RNA, and protein sequences. This tool provides comprehensive per-sequence statistics, formatted table output, and optional CSV export for bioinformatics workflows.

## Features

- **Multi-format Support**: Analyze DNA, RNA, and protein sequences
- **Comprehensive Analysis**: 
  - Sequence length calculation
  - GC content analysis (for DNA/RNA)
  - Amino acid composition (for proteins)
  - Invalid character detection
- **Flexible Input**: Process single files, multiple files, or entire directories
- **Recursive Processing**: Scan directories recursively for FASTA files
- **Multiple File Formats**: Supports `.fasta`, `.fa`, `.fna`, and `.faa` extensions
- **Formatted Output**: Pretty-printed table results in the terminal
- **CSV Export**: Optional CSV output for further analysis
- **Error Handling**: Robust validation and informative error messages

## Installation

### Requirements

- Python 3.7+
- No external dependencies (uses only Python standard library)

### Setup

1. Clone or download this repository
2. Ensure Python 3.7+ is installed on your system
3. Navigate to the project directory

## Usage

### Basic Usage

Analyze a single FASTA file:

```bash
python app.py data/genomic.fasta
```

### Analyze Multiple Files

```bash
python app.py data/*.fasta
```

### Specify Sequence Type

```bash
python app.py data/genomic.fasta --type DNA
```

**Supported types**: `DNA`, `RNA`, `PROTEIN`

### Recursive Directory Scanning

```bash
python app.py data/ --recursive
```

### Export Results to CSV

```bash
python app.py data/*.fasta --csv results.csv
```

### Combine Options

```bash
python app.py data/ --recursive --type DNA --csv analysis.csv
```

### View Help

```bash
python app.py --help
```

## File Structure

```
├── app.py                 # Main CLI entry point
├── fasta_analyzer.py      # Core analysis engine
├── data/                  # Sample FASTA files
│   ├── genomic.fasta
│   ├── mrna.fasta
│   └── protein.fasta
└── README.md              # This file
```

## Output Examples

### Terminal Output

```
File: genomic.fasta
┌────────────┬─────────┬──────────┬─────────────┐
│ Sequence   │ Length  │ Type     │ GC Content  │
├────────────┼─────────┼──────────┼─────────────┤
│ seq_001    │ 1,500   │ DNA      │ 45.2%       │
│ seq_002    │ 2,300   │ DNA      │ 52.1%       │
└────────────┴─────────┴──────────┴─────────────┘
```

### CSV Output

```
file_name,seq_id,seq_type,length,gc_content,invalid_chars
genomic.fasta,seq_001,DNA,1500,45.2,
genomic.fasta,seq_002,DNA,2300,52.1,
```

## Supported FASTA Formats

### DNA Sequences
Valid characters: A, T, G, C, N (ambiguous)

### RNA Sequences
Valid characters: A, U, G, C, N (ambiguous)

### Protein Sequences
Valid characters: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, X (ambiguous)

## Error Handling

The tool provides informative error messages for:
- **Missing files**: "File not found: [path]"
- **Empty files**: "File is empty: [path]"
- **Invalid FASTA format**: Missing header lines or malformed structure
- **Invalid characters**: Detected but reported with warnings
- **File permission errors**: Clear indication of access issues

## Examples

### Analyze Genomic Data

```bash
python app.py data/genomic.fasta --type DNA --csv genomic_analysis.csv
```

### Analyze All FASTA Files in Subdirectories

```bash
python app.py data/ --recursive --csv all_sequences.csv
```

### Analyze Protein Sequences with Default Type Detection

```bash
python app.py data/protein.fasta
```

## Technical Details

- **Language**: Python 3.7+
- **Dependencies**: None (standard library only)
- **Memory Efficient**: Streams large files without loading entirely into memory
- **Character Encoding**: UTF-8 with fallback error handling

## Use Cases

- **Genomics Research**: Analyze whole genomes and chromosome sequences
- **Metagenomics**: Process large collections of sequences from environmental samples
- **Protein Studies**: Analyze amino acid composition and properties
- **Data Quality**: Validate FASTA files and detect formatting issues
- **Batch Analysis**: Process multiple samples with consistent reporting

## License

This project is provided as-is for educational and research purposes.

## Contributing

Contributions are welcome! Feel free to:
- Report bugs or issues
- Suggest new features
- Improve documentation
- Optimize performance

## Support

For questions or issues, please refer to the error messages or check the file format of your FASTA files against the supported formats listed above.
