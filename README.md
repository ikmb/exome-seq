![](images/ikmb_bfx_logo.png)

# Exome-seq Pipeline

This pipeline offers a end-to-end workflow for exome analysis using several toolchains

- trimming with Fastp

- read alignment with BWA

- duplicate marking using Samtools

- vcf/gvcf calling with Deepvariant, Strelka and/or GATK

- Effect prediction with VEP and/or Haplosaurus

- SV calling using Manta

- CNV calling using CNVkit

- Expansion detection using ExpansionHunter

The result will be a multi-sample VCF file as well as a list of VCF files for each sample.

## Documentation 

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Benchmarking

Benchmarking of release 3.0  against genome-in-a-bottle NA12878

IDT xGEN, in-house:

| Variant  | Recall | Precision |
| -------- | ------ | --------- |
| Indel    | 0.92   | 0.98      |
| SNP      | 0.994  | 0.999     |

Agilent v7, external:

| Variant  | Recall | Precision |
| -------- | ------ | --------- |
| Indel    | 0.93   | 0.97      |
| SNP      | 0.99   | 0.998     |


