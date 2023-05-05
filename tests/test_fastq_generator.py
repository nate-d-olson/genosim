import pytest
from ..genosim.fastq_generator import Fastq, generate_read_sequences, generate_read_quality_string, generate_quality_strings, generate_read_names, generate_reads
from ..genosim.fasta_generator import Fasta


def test_fastq_init():
    fastq = Fastq()
    assert fastq.reads == []


def test_fastq_add_read():
    fastq = Fastq()
    fastq.add_read("test_read", "ATCG", "!#$%")
    assert fastq.reads == [("test_read", "ATCG", "!#$%")]


def test_generate_read_sequences():
    fasta = Fasta()
    fasta.add_sequence("test_sequence", "ATCGATCG")
    read_sequences, read_reference_names = generate_read_sequences(fasta, 5, 4, 0)
    assert len(read_sequences) == 5
    assert len(read_reference_names) == 5
    for seq in read_sequences:
        assert len(seq) == 4


def test_generate_read_quality_string():
    quality_string = generate_read_quality_string(10)
    assert len(quality_string) == 10


def test_generate_quality_strings():
    quality_strings = generate_quality_strings(5, 10, 0)
    assert len(quality_strings) == 5
    for qs in quality_strings:
        assert len(qs) == 10


def test_generate_read_names():
    read_names = generate_read_names("illumina", 5)
    assert len(read_names) == 5
    for rn in read_names:
        assert rn.startswith("@ILLUMINA")

    read_names = generate_read_names("pacbio", 5)
    assert len(read_names) == 5
    for rn in read_names:
        assert rn.startswith("@PACBIO_HIFI")

    read_names = generate_read_names("ont", 5)
    assert len(read_names) == 5
    for rn in read_names:
        assert rn.startswith("@ONT")

    with pytest.raises(ValueError):
        generate_read_names("invalid_technology", 5)


def test_generate_reads():
    fasta = Fasta()
    fasta.add_sequence("test_sequence", "ATCGATCG")

    reads_1, reads_2 = generate_reads(fasta, "illumina", 5, 4, 0)
    assert len(reads_1.reads) == 5
    assert len(reads_2.reads) == 5

    reads_1, reads_2 = generate_reads(fasta, "pacbio", 5, 4, 0)
    assert len(reads_1.reads) == 5
    assert reads_2 is None

    reads_1, reads_2 = generate_reads(fasta, "ont", 5, 4, 0)
    assert len(reads_1.reads) == 5
    assert reads_2 is None
