import random
import string
import os
import argparse
from typing import List, Tuple
from genosim.fasta_generator import Fasta

class Fastq:
    def __init__(self, reads=None):
        """
        Initialize the Fastq class with an optional list of reads.

        :param reads: A list of tuples containing read names, sequences, and quality scores. Default is None.
        :type reads: list[tuple[str, str, str]]
        """
        self.reads = reads if reads else []

    def add_read(self, read_name, sequence, quality_scores):
        """
        Add a read to the Fastq object.

        :param read_name: The read name.
        :type read_name: str
        :param sequence: The sequence.
        :type sequence: str
        :param quality_scores: The quality scores.
        :type quality_scores: str
        """
        self.reads.append((read_name, sequence, quality_scores))

    def parse_file(self, file_path):
        """
        Parse a Fastq file and load the reads into the Fastq object.

        :param file_path: The path of the Fastq file.
        :type file_path: str
        """
        with open(file_path, 'r') as file:
            lines = file.readlines()

            for i in range(0, len(lines), 4):
                read_name = lines[i].rstrip()[1:]
                sequence = lines[i + 1].rstrip()
                quality_scores = lines[i + 3].rstrip()
                self.reads.append((read_name, sequence, quality_scores))

    def write_file(self, file_path):
        """
        Write the reads in the Fastq object to a Fastq file.

        :param file_path: The path of the Fastq file to be written.
        :type file_path: str
        """
        with open(file_path, 'w') as file:
            for read_name, sequence, quality_scores in self.reads:
                file.write(f"@{read_name}\n")
                file.write(sequence + "\n")
                file.write("+\n")
                file.write(quality_scores + "\n")

def generate_read_sequences(reference_sequences: Fasta, num_reads: int, read_length: int, error_rate: float) -> List[str]:
    read_sequences = []
    read_reference_names = []

    for _ in range(num_reads):
        # Randomly select a reference sequence
        ref_name, ref_sequence = random.choice(reference_sequences.get_items())
        read_reference_names.append(ref_name)
        ref_length = len(ref_sequence)

        # Generate the start position for the read within the reference sequence
        start_pos = random.randint(0, ref_length - read_length)

        # Extract the read sequence from the reference sequence
        original_read = ref_sequence[start_pos:start_pos + read_length]

        # Introduce errors in the read sequence based on the error_rate
        read = [base if random.random() > error_rate else random.choice([b for b in 'ACGT' if b != base])
                for base in original_read]

        read_sequences.append(''.join(read))

    return read_sequences, read_reference_names


def generate_read_quality_string(read_length):
    return ''.join(random.choices(string.ascii_letters[26:60], k=read_length))

def generate_quality_strings(num_reads: int, read_length: int, error_rate: float) -> List[str]:
    quality_strings = []
    for _ in range(num_reads):
        quality_string = generate_read_quality_string(read_length)
        quality_strings.append(quality_string)
    return quality_strings

def generate_read_names(technology: str, num_reads: int) -> List[str]:
    read_names = []
    for i in range(num_reads):
        if technology.lower() == 'illumina':
            read_name = f"@ILLUMINA:test_read_{i}:8:1101:{i}:6789 1:N:0:ATCACG"
        elif technology.lower() == 'pacbio':
            read_name = f"@PACBIO_HIFI:test_read_{i}"
        elif technology.lower() == 'ont':
            read_name = f"@ONT:test_read_{i}"
        else:
            raise ValueError("Invalid technology specified. Choose from 'illumina', 'pacbio', or 'ont'.")
        read_names.append(read_name)
    return read_names

def generate_reads(reference_sequences: Fasta, technology: str, num_reads: int, read_length: int, error_rate: float) -> Tuple[Fastq, Fastq]:
    read_sequences, read_reference_names = generate_read_sequences(reference_sequences, num_reads, read_length, error_rate)
    quality_strings = generate_quality_strings(num_reads, read_length, error_rate)
    read_names = generate_read_names(technology, num_reads)

    if technology.lower() == 'illumina':
        reads_1 = Fastq()
        reads_2 = Fastq()

        for i in range(num_reads):
            reads_1.add_read(read_names[i], read_sequences[i], quality_strings[i])
            reads_2.add_read(read_names[i], read_sequences[i], quality_strings[i])

        return reads_1, reads_2

    else:
        reads = Fastq()

        for i in range(num_reads):
            reads.add_read(read_names[i], read_sequences[i], quality_strings[i])

        return reads, None

def generate_illumina_reads(reference_sequences: Fasta, num_pairs: int, read_length: int, insert_size: int, error_rate: float) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]], List[Tuple[str, str]]]:
    read_sequences_1 = generate_read_sequences(reference_sequences, num_pairs, read_length, error_rate)
    read_sequences_2 = [seq[-(insert_size - read_length):] + seq[:-(insert_size - read_length)] for seq in read_sequences_1]
    quality_strings = generate_quality_strings(num_pairs * 2, read_length, error_rate)
    quality_strings_1 = quality_strings[:num_pairs]
    quality_strings_2 = quality_strings[num_pairs:]
    read_names = generate_read_names('illumina', num_pairs * 2)
    read_names_1 = read_names[:num_pairs]
    read_names_2 = read_names[num_pairs:]
    paired_read_sequences = [(read1, read2) for read1, read2 in zip(read_sequences_1, read_sequences_2)]
    paired_quality_strings = [(qs1, qs2) for qs1, qs2 in zip(quality_strings_1, quality_strings_2)]
    paired_read_names = [(rn1, rn2) for rn1, rn2 in zip(read_names_1, read_names_2)]
    return paired_read_sequences, paired_quality_strings, paired_read_names


def generate_pacbio_reads(reference_sequences: Fasta, num_reads: int, read_length: int, error_rate: float) -> Tuple[List[str], List[str], List[str]]:
    return generate_read_sequences(reference_sequences, num_reads, read_length, error_rate)


def generate_ont_reads(reference_sequences: Fasta, num_reads: int, read_length: int, error_rate: float) -> Tuple[List[str], List[str], List[str]]:
    return generate_read_sequences(reference_sequences, num_reads, read_length, error_rate)

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate reads from a FASTA sequence.")
    parser.add_argument("-i", "--input", type=str, required=True, help="Input reference FASTA file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output FASTQ file.")
    parser.add_argument("-n", "--num_reads", type=int, default=1000, help="Number of reads to generate (default: 1000).")
    parser.add_argument("-l", "--read_length", type=int, default=150, help="Length of each read (default: 150).")
    parser.add_argument("-is", "--insert_size", type=int, default=300, help="Insert size for paired-end reads (default: 300).")
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed for reproducible results (default: None).")
    parser.add_argument("--technology", type=str, default="illumina", choices=["illumina", "pacbio", "ont"], help="Sequencing technology (default: 'illumina').")

    # Parse the arguments
    args = parser.parse_args()

    # Check if the output file already exists
    if os.path.exists(args.output):
        print(f"Error: Output file '{args.output}' already exists.")
        return

    # Set the random seed
    if args.seed is not None:
        random.seed(args.seed)

    # Read the input reference FASTA file
    ref = Fasta()
    ref.parse_file(args.input)

    # Read the input reference FASTA file
    ref = Fasta()
    ref.parse_file(args.input)

    # Select the function to generate reads based on the sequencing technology
    if args.technology == "illumina":
        error_rate = 0.001
    elif args.technology == "pacbio":
        error_rate = 0.01
    elif args.technology == "ont":
        error_rate = 0.1

    # Generate the reads
    reads_1, reads_2 = generate_reads(ref, args.technology, args.num_reads, args.read_length, error_rate)

    if reads_2 is None:  # Single-end reads
        reads_1.write_file(args.output)
        print(f"Generated {args.num_reads} reads from '{args.input}' and saved to '{args.output}'.")
    else:  # Paired-end reads
        output_prefix, output_ext = os.path.splitext(args.output)
        output_1 = f"{output_prefix}_1{output_ext}"
        output_2 = f"{output_prefix}_2{output_ext}"
        reads_1.write_file(output_1)
        reads_2.write_file(output_2)
        print(f"Generated {args.num_reads} paired-end reads from '{args.input}' and saved to '{output_1}' and '{output_2}'.")

if __name__ == "__main__":
    main()