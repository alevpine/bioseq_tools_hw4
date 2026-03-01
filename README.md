# Bioseq Tools

A Python toolkit for basic bioinformatics operations using OOP approach:
working with DNA, RNA, protein sequences, FASTQ filtering, and biological file processing.

## Overview

### Biological Sequence Classes

The core of the project is a hierarchy of classes for working with biological sequences:

- **BiologicalSequence** — abstract base class defining the interface: `len()`, indexing/slicing, `str()`, `check_alphabet()`.
- **NucleicAcidSequence** — base class for nucleic acids with `complement()`, `reverse()`, `reverse_complement()`.
  - **DNASequence** — DNA sequences, includes `transcribe()` method.
  - **RNASequence** — RNA sequences.
- **AminoAcidSequence** — protein sequences with `count_amino_acid()` method.

### FASTQ Filtering

Function `filter_fastq()` filters FASTQ files by GC content, read length, and average quality.
Implemented using **Biopython** (`SeqIO`, `gc_fraction`).


### Bio Files Processor

Additional tools for biological file formats:
- `convert_multiline_fasta_to_oneline()` — convert multiline FASTA to one-line format.
- `parse_blast_output()` — extract best hits from BLAST results.

## Installation

```bash
pip install -r requirements.txt
```

## Example Usage

### DNA / RNA

```python
from bioseq_tools import DNASequence, RNASequence

dna = DNASequence("ATGCGT")
print(dna)                      # ATGCGT
print(len(dna))                 # 6
print(dna[1:4])                 # TGC
print(dna.complement())         # TACGCA
print(dna.reverse_complement()) # ACGCAT
print(dna.transcribe())         # AUGCGU

rna = RNASequence("AUGCGU")
print(rna.complement())         # UACGCA
```

### Protein

```python
from bioseq_tools import AminoAcidSequence

protein = AminoAcidSequence("MVLSPADKTN")
print(len(protein))                  # 10
print(protein[0:3])                  # MVL
print(protein.count_amino_acid("V")) # 1
```

### FASTQ Filtering

```python
from bioseq_tools import filter_fastq

filter_fastq(
    input_fastq="example_data/sample.fastq",
    output_fastq="filtered_sample.fastq",
    gc_bounds=(40, 60),
    length_bounds=(50, 300),
    quality_threshold=30
)
```
### Bio Files Processor

```python
from bio_files_processor import convert_multiline_fasta_to_oneline, parse_blast_output

convert_multiline_fasta_to_oneline("input.fasta")
parse_blast_output("blast_results.txt", "parsed_hits.txt")
```

## Requirements

See `requirements.txt`.
