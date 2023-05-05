import pytest
import random
from io import StringIO
import os
from ..genosim.fasta_generator import Fasta
from ..genosim.fasta_generator import generate_reference_sequence
from ..genosim.fasta_generator import generate_multi_chromosome_reference


def test_fasta_init():
    fasta = Fasta()
    assert fasta.sequences == {}
    assert fasta.items == []

def test_fasta_add_sequence():
    fasta = Fasta()
    fasta.add_sequence("seq1", "ATGC")
    assert fasta.sequences == {"seq1": "ATGC"}
    assert fasta.items == [("seq1", "ATGC")]

def test_fasta_parse_file(tmp_path):
    file_path = tmp_path / "test_fasta.fasta"
    with open(file_path, "w") as file:
        file.write(">seq1\nATGC\n>seq2\nCGTA\n")
    fasta = Fasta()
    fasta.parse_file(file_path)
    assert fasta.sequences == {"seq1": "ATGC", "seq2": "CGTA"}
    assert fasta.items == [("seq1", "ATGC"), ("seq2", "CGTA")]

def test_fasta_write_file(tmp_path):
    file_path = tmp_path / "test_fasta.fasta"
    fasta = Fasta({"seq1": "ATGC", "seq2": "CGTA"})
    fasta.write_file(file_path)
    with open(file_path, "r") as file:
        content = file.read()
    assert content == ">seq1\nATGC\n>seq2\nCGTA\n"

def test_fasta_get_sequence():
    fasta = Fasta({"seq1": "ATGC", "seq2": "CGTA"})
    assert fasta.get_sequence("seq1") == "ATGC"
    assert fasta.get_sequence("seq2") == "CGTA"
    assert fasta.get_sequence("seq3") == None

def test_fasta_get_items():
    fasta = Fasta({"seq1": "ATGC", "seq2": "CGTA"})
    assert fasta.get_items() == [("seq1", "ATGC"), ("seq2", "CGTA")]


@pytest.mark.parametrize("length, homopolymer_length, repeat_length, repeat_count", [
    (20, 0, 0, 0),
    (20, 5, 0, 0),
    (20, 0, 3, 2),
    (20, 5, 3, 2),
])
def test_generate_reference_sequence(length, homopolymer_length, repeat_length, repeat_count):
    seq = generate_reference_sequence(length, homopolymer_length, repeat_length, repeat_count)
    assert len(seq) == length

def test_generate_multi_chromosome_reference():
    chromosomes = ["chr1", "chr2"]
    length = 10
    fasta = generate_multi_chromosome_reference(chromosomes, length, seed=42)

    assert len(fasta.sequences) == 2
    assert fasta.get_items()[0][0] == "chr1"
    assert len(fasta.get_sequence("chr1")) == length
    assert fasta.get_items()[1][0] == "chr2"
    assert len(fasta.get_sequence("chr2")) == length
