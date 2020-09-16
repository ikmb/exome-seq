![](images/ikmb_bfx_logo.png)

# Exome-seq Pipeline

This pipeline offers a end-to-end workflow for exome analysis using the DeepVariant toolchain

- trimming with Fastp

- read alignment with BWA

- duplicate marking using Samtools

- vcf/gvcf calling with Deepvariant

- joint variant calling with GLNexus

The result will be a multi-sample VCF file as well as a list of VCF files for each sample.

## Documentation 

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)
