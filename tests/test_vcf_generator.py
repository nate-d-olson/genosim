import os
import tempfile
import pysam
from ..genosim.vcf_generator import generate_vcf_from_bam

def test_generate_vcf_from_bam():
    # Use sample BAM and reference FASTA files for testing
    bam_file = "tests/data/test.bam"
    ref_fasta_file = "tests/data/test.fasta"

    # Create a temporary directory for the output files
    with tempfile.TemporaryDirectory() as tmpdir:
        output_vcf_file = os.path.join(tmpdir, "output.vcf.gz")

        # Run the generate_vcf_from_bam function
        generate_vcf_from_bam(bam_file, ref_fasta_file, output_vcf_file)

        # Check if the output VCF file and its index were generated
        assert os.path.isfile(output_vcf_file), "Output VCF file not generated"
        assert os.path.isfile(output_vcf_file + ".tbi"), "Output VCF index file not generated"

        # Check if the output VCF file can be opened with pysam
        try:
            with pysam.VariantFile(output_vcf_file) as vcf:
                assert len(vcf.header.samples) == 1, "Incorrect number of samples in VCF header"
                assert "sample1" in vcf.header.samples, "Sample1 not found in VCF header"
        except Exception as e:
            pytest.fail(f"Error opening output VCF file: {e}")
