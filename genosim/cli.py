import argparse
import os
import tempfile

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Generate a set of test files for a specified sequencing technology.")
    parser.add_argument("-t", "--technology", type=str, default="illumina", choices=["illumina", "pacbio", "ont"], help="Sequencing technology 
(default: 'illumina').")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output file prefix.")
    parser.add_argument("--num_chromosomes", type=int, default=2, help="Number of chromosomes to generate (default: 2).")
    parser.add_argument("--num_reads", type=int, default=1000, help="Number of reads to generate (default: 1000).")
    parser.add_argument("--read_length", type=int, default=150, help="Length of each read (default: 150).")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducible results (default: None).")

    # Parse the arguments
    args = parser.parse_args()

    # Generate a multi-chromosome FASTA file
    reference_file = args.output_prefix + "_reference.fasta"
    generate_multichromosome_fasta(args.num_chromosomes, reference_file, args.seed)

    # Generate reads
    reads_file = args.output_prefix + "_reads.fastq"
    generate_reads_cli = {
        "illumina": generate_illumina_reads,
        "pacbio": generate_pacbio_hifi_reads,
        "ont": generate_ont_reads
    }[args.technology]

    with tempfile.NamedTemporaryFile(mode="w") as temp_fasta:
        with open(reference_file, "r") as ref_fasta:
            temp_fasta.write(ref_fasta.read())
            temp_fasta.flush()

        generate_reads_from_fasta(temp_fasta.name, reads_file, args.num_reads, args.read_length, args.seed, generate_reads_cli)

    # Generate BAM file
    bam_file = args.output_prefix + "_alignment.bam"
    generate_bam_from_fastq(reads_file, reference_file, bam_file)

    # Generate VCF file
    vcf_file = args.output_prefix + "_variants.vcf"
    generate_vcf_from_bam(bam_file, reference_file, vcf_file)

    print(f"Generated test files:\n  Reference: {reference_file}\n  Reads: {reads_file}\n  Alignment: {bam_file}\n  Variants: {vcf_file}")

if __name__ == "__main__":
    main()

