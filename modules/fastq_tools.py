def is_nucleic_acid(sequence: str) -> bool:
    "Check DNA/RNA sequence."
    dna = {"A", "T", "G", "C"}
    rna = {"A", "U", "G", "C"}
    seq_upper = sequence.upper()
    return set(seq_upper).issubset(dna) or set(seq_upper).issubset(rna)


def average_quality(quality: str) -> float:
    "Calculate average quality (Phred33)."
    return sum(ord(char) - 33 for char in quality) / len(quality) if quality else 0


def gc_content(sequence: str) -> float:
    "Calculate GC content percentage."
    gc_count = sum(1 for nuc in sequence.upper() if nuc in "GC")
    return (gc_count / len(sequence)) * 100 if sequence else 0


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
) -> dict[str, tuple[str, str]]:
    "Filter FASTQ sequences by GC content, length, and average quality."

    def in_bounds(value: float, bounds: float | tuple[float, float]) -> bool:
        if isinstance(bounds, tuple):
            return bounds[0] <= value <= bounds[1]
        return value <= bounds

    filtered = {}
    for read_name, (sequence, quality) in seqs.items():
        if not is_nucleic_acid(sequence):
            print(f"Warning: Invalid sequence skipped: {sequence}")
            continue

        gc = gc_content(sequence)
        seq_length = len(sequence)
        avg_qual = average_quality(quality)

        if in_bounds(gc, gc_bounds) and in_bounds(seq_length, length_bounds) and avg_qual >= quality_threshold:
            filtered[read_name] = (sequence, quality)

    return filtered

