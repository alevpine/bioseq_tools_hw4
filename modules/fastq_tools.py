import os
from modules.dna_rna_tools import is_nucleic_acid


def average_quality(quality: str) -> float:
    """
    Calculate average quality score (Phred33).
    """
    if not quality:
        return 0
    return sum(ord(ch) - 33 for ch in quality) / len(quality)


def gc_content(sequence: str) -> float:
    """
    Calculate GC content percentage.
    """
    seq = sequence.lower()
    gc = seq.count("g") + seq.count("c")
    return (gc / len(seq)) * 100 if seq else 0


def in_bounds(value: float, bounds: float | tuple[float, float]) -> bool:
    """
    Check if a value is within given bounds.
    """
    if isinstance(bounds, tuple):
        return bounds[0] <= value <= bounds[1]
    return value <= bounds


def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: float | tuple[float, float] = (0, 100),
    length_bounds: int | tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
) -> dict[str, tuple[str, str]]:
    """
    Filter FASTQ sequences by GC content, length, and quality threshold.
    """
    filtered = {}

    for name, (seq, qual) in seqs.items():
        if not is_nucleic_acid(seq):
            continue

        gc = gc_content(seq)
        avg_q = average_quality(qual)
        length = len(seq)

        if (
            in_bounds(gc, gc_bounds)
            and in_bounds(length, length_bounds)
            and avg_q >= quality_threshold
        ):
            filtered[name] = (seq, qual)

    return filtered


def read_fastq(input_fastq: str) -> dict[str, tuple[str, str]]:
    """
    Read a FASTQ file and return a dictionary {read_name: (sequence, quality)}.
    """
    seqs = {}
    with open(input_fastq) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # skip "+"
            qual = f.readline().strip()
            seqs[header[1:]] = (seq, qual)
    return seqs


def write_fastq(seqs: dict[str, tuple[str, str]], output_fastq: str):
    """
    Write filtered sequences to a new FASTQ file in the "filtered" folder.
    """
    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq)

    base, ext = os.path.splitext(output_path)
    i = 1
    while os.path.exists(output_path):
        output_path = f"{base}_{i}{ext}"
        i += 1

    with open(output_path, "w") as out:
        for name, (seq, qual) in seqs.items():
            out.write(f"@{name}\n{seq}\n+\n{qual}\n")

