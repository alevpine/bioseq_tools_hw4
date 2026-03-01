from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# 1 Abstract Sequences


class BiologicalSequence(ABC):
    """Abstract base class for biological sequences."""

    @abstractmethod
    def __len__(self) -> int:
        """Return length of the sequence."""
        pass

    @abstractmethod
    def __getitem__(self, index):
        """Get element by index or slice."""
        pass

    @abstractmethod
    def __str__(self) -> str:
        """Return a nice string representation."""
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        """Check if the sequence alphabet is valid."""
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
        """Check that all characters belong to the valid alphabet."""
        if not self._alphabet:
            raise NotImplementedError(
                "Cannot check alphabet for base NucleicAcidSequence"
            )
        return set(self._sequence).issubset(self._alphabet)

    def complement(self):
        """Return complement sequence."""
        if not self._complement_map:
            raise NotImplementedError(
                "Cannot complement base NucleicAcidSequence"
            )
        comp = "".join(
            self._complement_map[nuc] for nuc in self._sequence
        )
        return type(self)(comp)

    def reverse(self):
        """Return reversed sequence."""
        return type(self)(self._sequence[::-1])

    def reverse_complement(self):
        """Return reverse complement sequence."""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences."""

    _alphabet = {"A", "T", "G", "C"}
    _complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        """Transcribe DNA to RNA. Returns RNASequence."""
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
        """Check that all characters are standard amino acids."""
        return set(self._sequence).issubset(self._alphabet)

    def count_amino_acid(self, amino_acid: str) -> int:
        """Count how many times a given amino acid appears."""
        return self._sequence.count(amino_acid.upper())


# 2 FastQ filter with Biopython


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

    Arguments:
        input_fastq: path to input FASTQ file.
        output_fastq: path to output FASTQ file.
        gc_bounds: tuple (min, max) for GC content in percent.
        length_bounds: tuple (min, max) for read length.
        quality_threshold: minimum average Phred33 quality score.
    """
    filtered_records = []

    for record in SeqIO.parse(input_fastq, "fastq"):
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
