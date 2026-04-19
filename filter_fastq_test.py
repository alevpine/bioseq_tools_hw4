import os
import pytest
import logging
from bioseq_tools import filter_fastq


def make_fastq_file(path: str, records: list[tuple[str, str, str]]):
    """
    Write records to a FASTQ file.
    Each record is (name, sequence, quality).
    """
    with open(path, "w") as f:
        for name, seq, qual in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


SAMPLE_RECORDS = [
    ("read1", "ATGCATGCAT", "IIIIIIIIII"),
    ("read2", "GCGCGCGCAT", "##########"),
    ("read3", "ATGC", "IIII"),
]


class TestFilterByGC:
    """Tests related to GC content filtering."""

    def test_gc_filter_keeps_matching(self, tmp_path):
        """Test that reads within GC bounds are kept."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file, gc_bounds=(40, 60))

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        names = [r.id for r in results]
        assert "read1" in names
        assert "read3" in names

    def test_gc_filter_removes_non_matching(self, tmp_path):
        """Test that reads outside GC bounds are removed."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file, gc_bounds=(40, 60))

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        names = [r.id for r in results]
        assert "read2" not in names


class TestFilterByLength:
    """Tests related to length filtering."""

    def test_length_filter_short_reads(self, tmp_path):
        """Test that short reads are filtered out by length bounds."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file, length_bounds=(5, 100))

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        names = [r.id for r in results]
        assert "read3" not in names
        assert len(results) == 2


class TestFilterByQuality:
    """Tests related to quality filtering."""

    def test_quality_filter_removes_low_quality(self, tmp_path):
        """Test that low quality reads are removed."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file, quality_threshold=30)

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        names = [r.id for r in results]
        assert "read2" not in names
        assert "read1" in names


class TestFilterFileIO:
    """Tests related to file reading and writing."""

    def test_output_file_is_created(self, tmp_path):
        """Test that output file is created after filtering."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file)

        assert os.path.exists(output_file)

    def test_all_reads_pass_with_no_filters(self, tmp_path):
        """Test that all reads pass when no filters are applied."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file)

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        assert len(results) == 3

    def test_empty_output_with_strict_filters(self, tmp_path):
        """Test that no reads pass with very strict filters."""
        input_file = str(tmp_path / "input.fastq")
        output_file = str(tmp_path / "output.fastq")
        make_fastq_file(input_file, SAMPLE_RECORDS)

        filter_fastq(input_file, output_file, quality_threshold=99)

        from Bio import SeqIO
        results = list(SeqIO.parse(output_file, "fastq"))
        assert len(results) == 0


class TestFilterErrors:
    """Tests related to error handling."""

    def test_missing_input_file_raises_error(self, tmp_path):
        """Test that FileNotFoundError is raised for missing input."""
        output_file = str(tmp_path / "output.fastq")

        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent_file.fastq", output_file)
