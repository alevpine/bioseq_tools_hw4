import argparse
import logging
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

logger = logging.getLogger("fastq_filter")


# Abstract Sequences


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences."""

    @abstractmethod
    def __len__(self) -> int:
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self) -> str:
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    """Base class for DNA and RNA sequences."""

    _alphabet: set = set()
    _complement_map: dict = {}

    def __init__(self, sequence: str):
        self._sequence = sequence.upper()
        if not self.check_alphabet():
            raise ValueError(
                f"Invalid characters in sequence for {type(self).__name__}"
            )

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index):
        result = self._sequence[index]
        if isinstance(index, slice):
            return type(self)(result)
        return result

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        return f"{type(self).__name__}('{self._sequence}')"

    def check_alphabet(self) -> bool:
        if not self._alphabet:
            raise NotImplementedError(
                "Cannot check alphabet for base NucleicAcidSequence"
            )
        return set(self._sequence).issubset(self._alphabet)

    def complement(self):
        if not self._complement_map:
            raise NotImplementedError(
                "Cannot complement base NucleicAcidSequence"
            )
        comp = "".join(
            self._complement_map[nuc] for nuc in self._sequence
        )
        return type(self)(comp)

    def reverse(self):
        return type(self)(self._sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences."""

    _alphabet = {"A", "T", "G", "C"}
    _complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        rna_seq = self._sequence.replace("T", "U")
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences."""

    _alphabet = {"A", "U", "G", "C"}
    _complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}


class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid (protein) sequences."""

    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence: str):
        self._sequence = sequence.upper()
        if not self.check_alphabet():
            raise ValueError("Invalid characters in protein sequence")

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index):
        result = self._sequence[index]
        if isinstance(index, slice):
            return AminoAcidSequence(result)
        return result

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        return f"AminoAcidSequence('{self._sequence}')"

    def check_alphabet(self) -> bool:
        return set(self._sequence).issubset(self._alphabet)

    def count_amino_acid(self, amino_acid: str) -> int:
        return self._sequence.count(amino_acid.upper())


# FastQ filter with Biopython


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: tuple[float, float] = (0, 100),
    length_bounds: tuple[int, int] = (0, 2**32),
    quality_threshold: float = 0,
) -> None:
    """
    Filter FASTQ file by GC content, length, and average quality.
    Uses Biopython for reading, writing, and calculations.
    """
    logger.info(
        f"Starting filtering: input={input_fastq}, output={output_fastq}, "
        f"gc_bounds={gc_bounds}, length_bounds={length_bounds}, "
        f"quality_threshold={quality_threshold}"
    )

    try:
        records = list(SeqIO.parse(input_fastq, "fastq"))
    except FileNotFoundError:
        logger.error(f"Input file not found: {input_fastq}")
        raise
    except Exception as error:
        logger.error(f"Error reading input file: {error}")
        raise

    filtered_records = []

    for record in records:
        seq_len = len(record.seq)
        gc = gc_fraction(record.seq) * 100
        qualities = record.letter_annotations["phred_quality"]
        avg_quality = sum(qualities) / len(qualities) if qualities else 0

        if (
            length_bounds[0] <= seq_len <= length_bounds[1]
            and gc_bounds[0] <= gc <= gc_bounds[1]
            and avg_quality >= quality_threshold
        ):
            filtered_records.append(record)

    SeqIO.write(filtered_records, output_fastq, "fastq")

    logger.info(
        f"Filtering done: {len(filtered_records)} of {len(records)} reads passed"
    )


# Command-line interface


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ file by GC content, length, and quality."
    )
    parser.add_argument("input", help="Path to input FASTQ file")
    parser.add_argument("output", help="Path to output FASTQ file")
    parser.add_argument(
        "--gc_bounds", nargs=2, type=float, default=[0, 100],
        help="Min and max GC content in percent (default: 0 100)"
    )
    parser.add_argument(
        "--length_bounds", nargs=2, type=int, default=[0, 2**32],
        help="Min and max read length (default: 0 4294967296)"
    )
    parser.add_argument(
        "--quality_threshold", type=float, default=0,
        help="Minimum average quality score (default: 0)"
    )

    args = parser.parse_args()

    logging.basicConfig(
        filename="fastq_filter.log",
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    filter_fastq(
        input_fastq=args.input,
        output_fastq=args.output,
        gc_bounds=tuple(args.gc_bounds),
        length_bounds=tuple(args.length_bounds),
        quality_threshold=args.quality_threshold,
    )


if __name__ == "__main__":
    main()
