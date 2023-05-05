import os
import tempfile
import pysam
from ..genosim.bam_generator import generate_bam_from_fastq

# Replace 'your_script_name' with the actual name of your script containing the generate_bam_from_fastq function

def test_generate_bam_from_fastq_single_end():
    # Setup
    fastq_file = "tests/data/test.fastq"
    ref_fasta_file = "tests/data/test.fasta"
    read_group_info = {"ID": "test", "SM": "toy"}
    with tempfile.NamedTemporaryFile(suffix=".bam") as temp_bam_file:
        output_bam_file = temp_bam_file.name

        # Call the function
        generate_bam_from_fastq([fastq_file], ref_fasta_file, output_bam_file, read_group_info, paired=False)

        # Check if the output BAM file is not empty
        with pysam.AlignmentFile(output_bam_file, "rb") as bam_file:
            assert sum(1 for _ in bam_file) > 0

        # Cleanup
        os.remove(output_bam_file.replace(".bam", "_sorted.bam"))
        os.remove(output_bam_file.replace(".bam", "_sorted.bam.bai"))

def test_generate_bam_from_fastq_paired_end():
    # Setup
    fastq_files = ["tests/data/test-ill_1.fastq", "tests/data/test-ill_2.fastq"]
    ref_fasta_file = "tests/data/test.fasta"
    read_group_info = {"ID": "test", "SM": "toy"}
    with tempfile.NamedTemporaryFile(suffix=".bam") as temp_bam_file:
        output_bam_file = temp_bam_file.name

        # Call the function
        generate_bam_from_fastq(fastq_files, ref_fasta_file, output_bam_file, read_group_info, paired=True)

        # Check if the output BAM file is not empty
        with pysam.AlignmentFile(output_bam_file, "rb") as bam_file:
            assert sum(1 for _ in bam_file) > 0

        # Cleanup
        os.remove(output_bam_file.replace(".bam", "_sorted.bam"))
        os.remove(output_bam_file.replace(".bam", "_sorted.bam.bai"))

