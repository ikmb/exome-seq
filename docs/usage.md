# Usage information

## Basic execution

The following command will execute this pipeline; the options will be discussed in the following.

If you are at the IKMB, the following will work:

`nextflow run ikmb/exome-seq --samples Samples.csv --assembly GRCh38 --kit xGen --email 'hello@gmail.com' 

If you try to run the pipeline on another system, you will need to configure a profile (see the installation instructions):

`nextflow run /path/to/main.nf --samples Samples.csv --assembly GRCh38 --kit xGen --email 'hello@gmail.com' -profile your_profile`


## Mandatory arguments

### `--samples`
Information to this pipeline is given in form of a CSV sample sheet. This is to allow relevant, and possibly quite important, metadata to be included (such
as meaningful sample names, sequencing center etc).

For convenience, we have included a simple script (bin/samplesheet_from_folder.rb) which accepts the path to a folder fill of PE exome data and automatically
writes a basic sample sheet. Obviously, it cannot derive meaningful sample and patient IDs from such information; this you would have to edit manually, if
you so choose.

```bash
ruby bin/samplesheet_from_folder.rb -f /path/to/foler > Samplesheet.csv`
```

### `--assembly` 
The following human genome assembly versions are supported on RZCluster (see resources.config on how this is set):

- GRCh37 (the 1000 genomes reference)
- GRCh38 (the current human reference assembly without ALT loci)
- hg19 (another version of GRCh37, also referred to as the UCSC reference)

### `--kit`
Each exome capture kit has a target and a bait definition, i.e. information about the exons it enriches and the specific RNA bait sequences that are used 
for capture. This information is important so the pipeline knows whichs regions of the genome to analyze and how to compute run metrics. 

We have included these files for the following kits and genome assemblies:

`Nextera` [hg19, GRCH37, GRCh38]

`xGen` (original IDT xGen kit) [hg19, GRCh37, GRCh38]

`xGen_v2` (v2 release of the IDT xGen kit) [GRCh38]

`Pan_cancer` [hg19]

### `--email`
Your Email address in quotes to which the pipeline report is sent upon completion. 

## Optional arguments

## `--baits` | `--targets` 
If you have used any other type of kit for your enrichment, you are able to provide these from the command line during execution using `--baits` and 
`--targets`, respectively. Please not that these files must be in the Picard 
[interval_list](https://gatkforums.broadinstitute.org/gatk/discussion/1319/collected-faqs-about-interval-lists) format and have to be matched 
to the genome assembly (i.e. must have identical dictionary headers). 

## `--panel`
For practical reasons, it can be desirable to determine the coverage of a discrete set of target genes, such as for a gene panel. The pipeline currently 
supports the following panels:

- Dilatative Kardiomyopathie [cardio_dilatative]
- Hypertrophe Kardiomyopathie [cardio_hypertrophic]
- Non-Compaction Kardiomyopathie [cardio_non_compaction]
- Immune defect AK deficency (25kb panel) [ gene_immunedefect_ak-25kb ]
- Immune defect CVID (25kb panel) [ gene_immunedefect_cvid-25kb ]
- Immune defect full [ gene_immunedefect_full ]
- Immune defect with intestinal component (25kb panel) [ gene_immunedefect_with_intestine-25kb ]
- Breast cancer panel [ breast_cancer ]

Please not that this will also create additional run metrics, including a per-sample list of target exons that fall below a minimum sequence coverage. 

## `--max_length` 
Set this to a positive number to trim all reads down to a desired size. Default: no size-trimming.

## `--interval_padding`
Set this to a positive number to include flanking regions of exon targets in the analysis. Default: 10

## `--skip_multiqc`
Skip the sending of a QC report. Default: false

## `--cram`
Create CRAM instead of BAM files to save space. Not that CRAM files are slower to read by IGV. 

## `--vqsr`
Perform the variant score recalibration filtering workflow. This requires > 30 exomes to be analysed in parallel and is deactivated by default. 

## `--run_name`
Give this run a meaningful name (like a LIMS or project ID)

## `--fasta`
Provide path to a genome sequence in FASTA format (default: false, uses a pre-configured genome)

## `--dict`
Provide path to a genome sequence dictionary file (default: false, uses a pre-configured dictionary)

## `--dbsnp`
Provide path to a dbSNP reference VCF file for variant filtering (default: false, uses a pre-configured reference)

## `--g1k`
Provide path to a 1000genomes reference VCF file for variant filtering (default: false, uses a pre-configured reference)

## `--mills_indels`
Provide path to a indel reference VCF file for variant filtering (default: false, uses a pre-configured reference)

## `--omni`
Provide path to a SNP reference VCF file for variant filtering (default: false, uses a pre-configured reference)

## `--hapmap`
Provide path to a SNP reference VCF file for variant filtering (default: false, uses a pre-configured reference)

