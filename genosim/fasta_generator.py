import random
import os
import argparse

class Fasta:
    def __init__(self, sequences=None):
        """
        Initialize the Fasta class with an optional dictionary of sequences.

        :param sequences: A dictionary containing sequence names as keys and sequences as values. Default is None.
        :type sequences: dict[str, str]
        """
        self.sequences = sequences if sequences else {}
        self.items = list(self.sequences.items())

    def add_sequence(self, sequence_name, sequence):
        """
        Add a sequence to the Fasta object.

        :param sequence_name: The sequence name.
        :type sequence_name: str
        :param sequence: The sequence.
        :type sequence: str
        """
        self.sequences[sequence_name] = sequence
        self.items.append((sequence_name, sequence))

    def parse_file(self, file_path):
        """
        Parse a Fasta file and load the sequences into the Fasta object.

        :param file_path: The path of the Fasta file.
        :type file_path: str
        """
        with open(file_path, 'r') as file:
            sequence_name = None
            sequence = []

            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if sequence_name:
                        self.add_sequence(sequence_name, ''.join(sequence))
                    sequence_name = line[1:]
                    sequence = []
                else:
                    sequence.append(line)

            if sequence_name:
                self.add_sequence(sequence_name, ''.join(sequence))

    def write_file(self, file_path):
        """
        Write the sequences in the Fasta object to a Fasta file.

        :param file_path: The path of the Fasta file to be written.
        :type file_path: str
        """
        with open(file_path, 'w') as file:
            for sequence_name, sequence in self.items:
                file.write(f">{sequence_name}\n")
                file.write(sequence + "\n")

    def get_sequence(self, sequence_name):
        """
        Get a sequence from the Fasta object by its name.

        :param sequence_name: The sequence name.
        :type sequence_name: str
        :return: The sequence corresponding to the given name.
        :rtype: str
        """
        return self.sequences.get(sequence_name, None)

    def get_items(self):
        """
        Get the list of items (sequence name, sequence) in the Fasta object.

        :return: The list of items.
        :rtype: list[tuple[str, str]]
        """
        return self.items


def generate_reference_sequence(length, seed, homopolymer_length=0, repeat_length=0, repeat_count=0):
    """
    Generates a reference sequence of the desired length, with options to include specific sequence contexts and edge cases.

    :param length: The desired length of the reference sequence.
    :type length: int
    :param homopolymer_length: The length of a homopolymer to be included in the sequence. Default is 0 (no homopolymer).
    :type homopolymer_length: int
    :param repeat_length: The length of a repeat to be included in the sequence. Default is 0 (no repeats).
    :type repeat_length: int
    :param repeat_count: The number of times the repeat should occur in the sequence. Default is 0 (no repeats).
    :type repeat_count: int
    :return: A reference sequence of the specified length with the specified sequence contexts and edge cases.
    :rtype: str
    """
    nucleotides = ['A', 'C', 'G', 'T']
    random.seed = seed
    # Create a homopolymer if specified
    homopolymer = ""
    if homopolymer_length > 0:
        base = random.choice(nucleotides)
        homopolymer = base * homopolymer_length
    
    # Create a repeat sequence if specified
    repeat = ""
    if repeat_length > 0 and repeat_count > 0:
        repeat = ''.join(random.choices(nucleotides, k=repeat_length)) * repeat_count
    
    # Generate the random sequence
    random_length = length - (homopolymer_length + repeat_length * repeat_count)
    random_sequence = ''.join(random.choices(nucleotides, k=random_length))
    
    # Combine the sequence components
    sequence = random_sequence + homopolymer + repeat
    sequence = ''.join(random.sample(sequence, len(sequence)))

    return sequence

def generate_multi_chromosome_reference(chromosomes, length, seed=None):
    """
    Generate a multi-chromosome reference sequence.

    Args:
        chromosomes (List[str]): List of chromosome names.
        length (int): Length of each chromosome.
        seed (int, optional): Random seed for reproducible results. Defaults to None.

    Returns:
        Fasta: Fasta object representing the multi-chromosome reference sequence.
    """
    fasta = Fasta()
    for chromosome in chromosomes:
        sequence = generate_reference_sequence(length, seed)
        fasta.add_sequence(chromosome, sequence)
    return fasta


def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate a multi-chromosome reference sequence in FASTA format.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output FASTA file name.")
    parser.add_argument("-c", "--chromosomes", type=int, default=2, help="Number of chromosomes (default: 2).")
    parser.add_argument("-l", "--length", type=int, default=1000000, help="Length of each chromosome (default: 1000000).")
    parser.add_argument("-s", "--seed", type=int, default=None, help="Random seed for reproducible results (default: None).")

    # Parse the arguments
    args = parser.parse_args()

    # Check if the output file already exists
    if os.path.exists(args.output):
        print(f"Error: Output file '{args.output}' already exists.")
        return

    # Generate the multi-chromosome reference sequence
    chrom_names = [ f"chr{i}" for i in range(args.chromosomes)]
    reference_sequence = generate_multi_chromosome_reference(chrom_names, args.length, args.seed)

    # Write the reference sequence to the output FASTA file
    reference_sequence.write_file(args.output)

    print(f"Generated multi-chromosome reference sequence saved to '{args.output}'.")

if __name__ == "__main__":
    main()