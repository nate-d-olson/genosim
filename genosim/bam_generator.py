import argparse
import pysam
import array
from tqdm.auto import tqdm

def generate_bam_from_fastq(fastq_files, ref_fasta_file, output_bam_file, read_group_info, paired=False):
    """
    Generate a BAM file from FASTQ files and a reference FASTA file.

    :param fastq_files: A list of paths to the input FASTQ files containing reads.
    :type fastq_files: list[str]
    :param ref_fasta_file: The path to the reference FASTA file.
    :type ref_fasta_file: str
    :param output_bam_file: The path to the output BAM file.
    :type output_bam_file: str
    :param read_group_info: The read group information.
    :type read_group_info: dict
    :param paired: A flag indicating whether the input reads are paired-end. Default is False.
    :type paired: bool
    """
    # Load the reference FASTA file
    ref_fasta = pysam.FastaFile(ref_fasta_file)

    # Create a header for the SAM file
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [{"SN": ref, "LN": length} for ref, length in zip(ref_fasta.references, ref_fasta.lengths)],
        "RG": [read_group_info],
    }

    # Create an empty SAM file with the header information
    with pysam.AlignmentFile(output_bam_file, "wb", header=header) as output_sam:
        # Parse the input FASTQ files
        input_fastqs = [pysam.FastxFile(f) for f in fastq_files]
        
        # Create a progress bar
        progress_bar = tqdm(desc="Processing reads", unit="reads")
        while True:
            try:
                entries = [next(f) for f in input_fastqs]

                for i, entry in enumerate(entries):
                    # Create an alignment object for each read
                    aln = pysam.AlignedSegment()
                    aln.query_name = entry.name
                    aln.query_sequence = entry.sequence
                    aln.flag = 77 if paired and i == 0 else 141 if paired and i == 1 else 4  # Paired and unmapped reads
                    aln.reference_id = -1
                    aln.reference_start = -1
                    aln.mapping_quality = 0
                    aln.cigar = None
                    aln.next_reference_id = -1
                    aln.next_reference_start = -1
                    aln.template_length = 0
                    aln.query_qualities = array.array("B", (ord(c) - 33 for c in entry.quality))
                    aln.set_tag("RG", read_group_info["ID"])

                    # Add the alignment to the SAM file
                    output_sam.write(aln)

                # Update the progress bar
                progress_bar.update(len(entries))

            except StopIteration:
                break

        # Close the progress bar
        progress_bar.close()

        # Close the input FASTQ files
        for f in input_fastqs:
            f.close()

    # Sort and index the BAM file
    pysam.sort("-o", output_bam_file.replace(".bam", "_sorted.bam"), output_bam_file)
    pysam.index(output_bam_file.replace(".bam", "_sorted.bam"))

def main():
    parser = argparse.ArgumentParser(description="Generate a BAM file from FASTQ files and a reference FASTA file.")
    parser.add_argument("-f", "--fastq_files", nargs="+", required=True, type=argparse.FileType('r'), help="Input FASTQ files containing reads.")
    parser.add_argument("-r", "--ref_fasta_file", required=True, type=argparse.FileType('r'), help="Reference FASTA file.")
    parser.add_argument("-o", "--output_bam_file", required=True, help="Output BAM file (with .bam extension).")
    parser.add_argument("-g", "--read_group_info", default="ID:test,SM:toy", help="Read group information in the format 'ID:value,SM:value'. Default is 'ID:test,SM:toy'.")
    parser.add_argument("--paired", action="store_true", help="Indicate that the input reads are paired-end.")

    args = parser.parse_args()

    # Check output file extension
    if not args.output_bam_file.endswith(".bam"):
        parser.error("Output BAM file must have a .bam extension.")

    # Validate read group information
    try:
        read_group_info = {key_value.split(':')[0]: key_value.split(':')[1] for key_value in args.read_group_info.split(',')}
        if 'ID' not in read_group_info or 'SM' not in read_group_info:
            raise ValueError
    except ValueError:
        parser.error("Invalid read group information. Please use the format 'ID:value,SM:value'.")

    # Validate the number of input FASTQ files
    if args.paired and len(args.fastq_files) != 2:
        parser.error("Paired-end mode requires exactly 2 input FASTQ files.")
    elif not args.paired and len(args.fastq_files) != 1:
        parser.error("Single-end mode requires exactly 1 input FASTQ file.")

    # Close argparse.FileType instances to avoid file locks
    for f in args.fastq_files:
        f.close()
    args.ref_fasta_file.close()

    # Get the file paths from argparse.FileType instances
    fastq_files = [f.name for f in args.fastq_files]
    ref_fasta_file = args.ref_fasta_file.name

    print("Processing input files...")
    generate_bam_from_fastq(fastq_files, ref_fasta_file, args.output_bam_file, read_group_info, args.paired)
    print("Done!")

if __name__ == "__main__":
    main()
