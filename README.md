# Bioseq Tools (HW5 – Files & Modules)

A simple Python toolkit for basic bioinformatics operations:  
working with DNA/RNA sequences, FASTQ files, and other biological data formats.

---

## Overview

### DNA/RNA Tools
Functions for working with nucleotide sequences:
- `is_nucleic_acid(sequence: str) -> bool` — check if a sequence is valid DNA or RNA  
- `transcribe(sequence: str) -> str` — convert DNA to RNA  
- `reverse(sequence: str) -> str` — reverse sequence order  
- `complement(sequence: str) -> str` — get complementary sequence  
- `reverse_complement(sequence: str) -> str` — get reverse complement  
- `run_dna_rna_tools(*args: str)` — main wrapper for all DNA/RNA operations  

---

### FASTQ Tools
Functions for reading, filtering, and writing FASTQ files:
- `read_fastq(input_fastq: str)` — read FASTQ file into a dictionary  
- `filter_fastq(seqs: dict, gc_bounds, length_bounds, quality_threshold)` — filter reads by GC content, length, and average quality  
- `write_fastq(seqs: dict, output_fastq: str)` — save filtered reads into `filtered/` folder, safely avoiding overwriting  
- `run_fastq_filter(input_fastq: str, output_fastq: str, **kwargs)` — main wrapper combining all above steps  

---

### Main Script — bioseq_tools.py
The main interface that connects both tools:
- imports and exposes DNA/RNA and FASTQ functions  
- allows you to run all operations from a single entry point  

Main functions:
- `run_dna_rna_tools(*args)` — for sequence analysis and transformations  
- `run_fastq_filter(input_fastq, output_fastq, **kwargs)` — for FASTQ file filtering and saving  

---

### Bio Files Processor
Additional tools for other biological formats:
- `convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str | None)` — make FASTA sequences single-line  
- `parse_blast_output(input_file: str, output_file: str)` — extract best hits from BLAST text results and sort them alphabetically  

## Example Usage

### DNA/RNA Tools

```python
from bioseq_tools import run_dna_rna_tools

# Check if sequence is valid DNA/RNA
print(run_dna_rna_tools("ATGC", "is_nucleic_acid"))      # True

# Transcribe DNA → RNA
print(run_dna_rna_tools("ATGC", "transcribe"))           # 'AUGC'

# Reverse, complement, and reverse complement
print(run_dna_rna_tools("ATGC", "reverse"))              # 'CGTA'
print(run_dna_rna_tools("ATGC", "complement"))           # 'TACG'
print(run_dna_rna_tools("ATGC", "reverse_complement"))   # 'GCAT'

# Multiple sequences
print(run_dna_rna_tools("ATG", "AAT", "reverse"))        # ['GTA', 'TAA']
````

After running these commands,
you will see the results of each DNA/RNA operation printed directly in the console.

---

### FASTQ Tools

```python
from bioseq_tools import run_fastq_filter

# Read, filter, and save FASTQ data automatically
run_fastq_filter(
    input_fastq="example_data/sample.fastq",
    output_fastq="filtered_sample.fastq",
    gc_bounds=(40, 60),
    quality_threshold=30
)
# Output file: filtered/filtered_sample.fastq
```

This will:

1. Read the FASTQ file,
2. Filter sequences by GC content and average quality,
3. Save the result safely to the `filtered/` folder.

---

### Main Script — bioseq_tools.py

```python
from bioseq_tools import run_dna_rna_tools, run_fastq_filter

# DNA/RNA operations
print(run_dna_rna_tools("ATGC", "reverse_complement"))   # 'GCAT'

# FASTQ filtering
run_fastq_filter(
    input_fastq="example_data/sample.fastq",
    output_fastq="filtered_output.fastq",
    gc_bounds=(40, 70),
    quality_threshold=25
)
```

`bioseq_tools.py` acts as a single entry point —
it lets you use both DNA/RNA and FASTQ functions from one file.

---

### Bio Files Processor

```python
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output

# Convert multiline FASTA → one-line format
convert_multiline_fasta_to_oneline("example_data/input.fasta")
# Creates: example_data/input_oneline.fasta

# Parse BLAST results and extract best hits
parse_blast_output("example_data/blast_results.txt", "parsed_hits.txt")
# Creates: parsed_hits.txt (sorted alphabetically)
```

These utilities handle biological text files directly —
you can quickly clean FASTA sequences or extract key results from BLAST outputs.

---

### Full Workflow Example

```python
from bioseq_tools import run_dna_rna_tools, run_fastq_filter
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output

# 1. DNA/RNA operations
dna = "ATGCGT"
print("Original:", dna)
print("Transcribed:", run_dna_rna_tools(dna, "transcribe"))
print("Reverse complement:", run_dna_rna_tools(dna, "reverse_complement"))

# 2. FASTQ filtering
run_fastq_filter(
    input_fastq="example_data/sample.fastq",
    output_fastq="filtered_sample.fastq",
    gc_bounds=(40, 60),
    quality_threshold=30
)

# 3. FASTA and BLAST file processing
convert_multiline_fasta_to_oneline("example_data/input.fasta")
parse_blast_output("example_data/blast_results.txt", "parsed_hits.txt")

print("All processing steps completed successfully.")
```

This single script demonstrates the entire workflow:

1. Works with DNA/RNA sequences,
2. Filters FASTQ reads,
3. Processes FASTA and BLAST files.
   All output files are saved automatically inside the project folders.

```

