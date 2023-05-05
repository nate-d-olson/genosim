import argparse
import os
import sys
import pysam

def generate_vcf_from_bam(bam_file, ref_fasta_file, output_vcf_file):
    """
    Generate a compressed VCF file from a BAM file and a reference FASTA file.

    :param bam_file: The path to the input BAM file.
    :type bam_file: str
    :param ref_fasta_file: The path to the reference FASTA file.
    :type ref_fasta_file: str
    :param output_vcf_file: The path to the output compressed VCF file.
    :type output_vcf_file: str
    """
    # Read the reference FASTA file
    ref_fasta = pysam.FastaFile(ref_fasta_file)

    # Create a header for the VCF file
    header = pysam.VariantHeader()
    for ref, length in zip(ref_fasta.references, ref_fasta.lengths):
        header.add_line(f'##contig=<ID={ref},length={length}>')

    header.add_sample("sample1")
    header.add_meta("FORMAT", items=[("ID", "GT"), ("Number", "1"), ("Type", "String"), ("Description", "Genotype")])
    header.add_meta("FORMAT", items=[("ID", "DP"), ("Number", "1"), ("Type", "Integer"), ("Description", "Read depth")])

    # Open the input BAM file
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Create an output compressed VCF file with the header information
        with pysam.VariantFile(output_vcf_file, "w", header=header) as output_vcf:
            # Iterate over the reference sequences in the BAM file
            for ref, length in zip(bam.references, bam.lengths):
                # Fetch the reads aligned to the current reference sequence
                for pileupcolumn in bam.pileup(ref):
                    # Perform variant calling (e.g., based on a simple threshold)
                    if pileupcolumn.nsegments > 0 and pileupcolumn.nsegments / length > 0.5:
                        var = output_vcf.new_record(contig=ref, pos=pileupcolumn.pos + 1)
                        var.ref = ref_fasta.fetch(ref, pileupcolumn.pos, pileupcolumn.pos + 1)
                        var.alts = ["<NON_REF>"]
                        var.samples["sample1"]["GT"] = (0, 1)
                        var.samples["sample1"]["DP"] = pileupcolumn.nsegments

                        # Add the variant to the output VCF file
                        output_vcf.write(var)

    # Index the output compressed VCF file
    pysam.tabix_index(output_vcf_file, preset='vcf')



def main():
    parser = argparse.ArgumentParser(description='Generate a compressed VCF file from a BAM file and a reference FASTA file.')
    parser.add_argument('-b', '--bam_file', required=True, help='Path to the input BAM file.')
    parser.add_argument('-r', '--ref_fasta_file', required=True, help='Path to the reference FASTA file.')
    parser.add_argument('-o', '--output_vcf_file', required=True, help='Path to the output compressed VCF file.')

    args = parser.parse_args()

    if not os.path.isfile(args.bam_file):
        print(f"Error: BAM file '{args.bam_file}' not found.")
        sys.exit(1)

    if not os.path.isfile(args.ref_fasta_file):
        print(f"Error: Reference FASTA file '{args.ref_fasta_file}' not found.")
        sys.exit(1)

    output_dir = os.path.dirname(args.output_vcf_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    generate_vcf_from_bam(args.bam_file, args.ref_fasta_file, args.output_vcf_file)
    print(f"Compressed VCF file generated at '{args.output_vcf_file}'.")

if __name__ == '__main__':
    main()
