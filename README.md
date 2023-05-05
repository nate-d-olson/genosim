# Bioinformatics Test Data Generator

The **Bioinformatics Test Data Generator** is a Python package that provides utilities for generating test datasets for use in developing 
bioinformatic pipelines for analyzing human genome sequencing data. It allows you to generate multi-chromosome reference sequences, simulated 
sequencing reads, alignments (BAM files), and variant calls (VCF files) based on different sequencing technologies.

## Features

- Generate multi-chromosome reference sequences in FASTA format.
- Simulate reads based on different sequencing technologies, such as Illumina, PacBio, and Oxford Nanopore Technologies.
- Generate alignments (BAM files) from simulated reads and reference sequences.
- Generate variant calls (VCF files) from alignments.

## Installation

You can install the package using pip:

```shell
pip install genosim
```


## Usage

The package provides a command-line interface (CLI) for generating the test datasets. Here are some example commands:

Generate a multi-chromosome reference sequence:

```shell
genosim --technology illumina --output_prefix my_test_data --num_chromosomes 2
```


Generate sequencing reads from a reference FASTA:

```shell
genosim --technology pacbio --output_prefix my_test_data --num_reads 1000 --read_length 100
```

Generate a BAM file from reads and reference:
```shell
genosim --technology ont --output_prefix my_test_data --read_group "RG1:ONT:1:lib1:sample1"
```


For more information and options, use the `--help` flag with the respective command.

## Documentation

For detailed usage instructions and API documentation, please refer to the 
[documentation](https://github.com/yourusername/bioinformatics_test_data_generator).

## Contributing

Contributions are welcome! Please see the [contribution guidelines](CONTRIBUTING.md) for more information.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.


