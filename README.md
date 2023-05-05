+![](images/ikmb_bfx_logo.png)

# Exome-seq Pipeline

This pipeline offers a end-to-end workflow for exome analysis using several toolchains

- trimming with [fastp](https://github.com/OpenGene/fastp)

- short read alignment using [BWA](https://github.com/lh3/bwa), [BWA2](https://github.com/bwa-mem2/bwa-mem2) or [Dragmap](https://github.com/Illumina/DRAGMAP)

- duplicate marking using [Samtools](https://github.com/samtools/samtools)

- germline SNP/INDEL calling with [Deepvariant](https://github.com/google/deepvariant), [Strelka](https://github.com/Illumina/strelka) and/or [GATK](https://github.com/broadinstitute/gatk)

- somatic SNP/INDEL calling with [Mutect2](https://github.com/broadinstitute/gatk) and [Strelka](https://github.com/Illumina/strelka)

- variant effect prediction with [VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and/or [Haplosaurus](https://www.ensembl.org/info/docs/tools/vep/haplo/index.html)

- protein-level effect prediction with [Haplosaurus](https://www.ensembl.org/info/docs/tools/vep/haplo/index.html) and [BCFtools](https://samtools.github.io/bcftools/howtos/csq-calling.html)

- germline and somatic SV calling using [Manta](https://github.com/Illumina/manta)

- germline and somatic CNV calling using [CNVkit](https://github.com/etal/cnvkit)

- Repeat expansion detection using [ExpansionHunter](https://github.com/Illumina/ExpansionHunter)

The result will be a multi-sample VCF file as well as a list of VCF files for each sample.

## Documentation 

1. [What happens in this pipeline?](docs/pipeline.md)
2. [Installation and configuration](docs/installation.md)
3. [Running the pipeline](docs/usage.md)
4. [Output](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Benchmarking

Benchmarking of release 4.3  against genome-in-a-bottle NA12878

### DeepVariant

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


