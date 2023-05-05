import random
import string
from collections import namedtuple

def generate_test_dataset(read_type, read_length, error_rate, num_reads, reference_sequence):
    """
    Generate a test dataset with specified characteristics.

    :param read_type: The read type, either 'single' or 'paired'.
    :type read_type: str
    :param read_length: The read length.
    :type read_length: int
    :param error_rate: The base error rate.
    :type error_rate: float
    :param num_reads: The number of reads to generate.
    :type num_reads: int
    :param reference_sequence: The reference sequence as a string.
    :type reference_sequence: str
    :return: A list of generated reads.
    :rtype: list[tuple[str, str]] for single-end reads or list[tuple[str, str, str, str]] for paired-end reads
    """
    if read_type not in ['single', 'paired']:
        raise ValueError("Invalid read_type. Must be 'single' or 'paired'.")

    Read = namedtuple("Read", "read_name seq1 seq2 qual1 qual2")

    def generate_read(read_length, error_rate, ref_seq):
        mutated_seq = []
        for base in ref_seq:
            if random.random() < error_rate:
                mutated_base = random.choice([b for b in 'ATCG' if b != base])
            else:
                mutated_base = base
            mutated_seq.append(mutated_base)
        return ''.join(mutated_seq)

    def generate_quality_string(read_length):
        return ''.join(random.choices(string.ascii_letters[26:60], k=read_length))

    generated_reads = []

    for i in range(num_reads):
        read_name = f"test_read_{i}"
        start_pos = random.randint(0, len(reference_sequence) - read_length)
        ref_seq = reference_sequence[start_pos: start_pos + read_length]
        seq1 = generate_read(read_length, error_rate, ref_seq)
        qual1 = generate_quality_string(read_length)

        if read_type == 'single':
            generated_reads.append(Read(read_name, seq1, None, qual1, None))
        else:
            seq2 = generate_read(read_length, error_rate, ref_seq)
            qual2 = generate_quality_string(read_length)
            generated_reads.append(Read(read_name, seq1, seq2, qual1, qual2))

    return generated_reads

