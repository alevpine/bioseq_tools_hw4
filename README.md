# Bioseq Tools (Python HW4, Modules)

A simple command-line utility for working with DNA/RNA sequences and FASTQ-like data.

## Features

### DNA/RNA Tools
- Check if a sequence is valid DNA or RNA (`is_nucleic_acid`)
- Transcribe DNA to RNA (`transcribe`)
- Reverse a sequence (`reverse`)
- Get complementary sequence (`complement`)
- Get reverse complementary sequence (`reverse_complement`)
- Main function `run_dna_rna_tools` works with one or multiple sequences

### FASTQ Tools
- Filter sequences by GC content, length, and average quality (`filter_fastq`)
- Supports Phred33 quality scores
- Returns only sequences that meet all specified criteria

## How It Works

### DNA/RNA Tools
`run_dna_rna_tools` accepts one or more sequences and the procedure name (last argument).  
For each sequence:
- Validity is checked (if procedure is not `is_nucleic_acid`)
- The corresponding procedure is executed
- Returns a single result for one sequence or a list for multiple sequences

### FASTQ Tools
`filter_fastq` accepts a dictionary of sequences:
- Key: sequence name  
- Value: tuple of `(sequence, quality)`  
Filters sequences according to GC content, length, and average quality thresholds.

## Example Usage

```python
from bioseq_tools import run_dna_rna_tools, filter_fastq

# DNA/RNA Tools
print(run_dna_rna_tools("TTUU", "is_nucleic_acid"))         # False
print(run_dna_rna_tools("ATG", "transcribe"))              # 'AUG'
print(run_dna_rna_tools("ATG", "reverse"))                 # 'GTA'
print(run_dna_rna_tools("AtG", "complement"))              # 'TaC'
print(run_dna_rna_tools("ATg", "reverse_complement"))      # 'cAT'
print(run_dna_rna_tools("ATG", "aT", "reverse"))           # ['GTA', 'Ta']

# FASTQ Tools
seqs = {
    "read1": ("ATGC", "IIII"),
    "read2": ("ATUU", "####"),
}
filtered = filter_fastq(seqs, gc_bounds=(50,100), quality_threshold=20)
print(filtered)  # Only sequences passing the filters
