# Usage information

## Basic execution

The following command will execute this pipeline; the options will be discussed in the following.

If you are at the IKMB, the following will work:

`nextflow run ikmb/exome-seq --samples Samples.csv --assembly GRCh38 --kit xGen --email 'hello@gmail.com`

If you need to run a specific "release" of the pipeline, you can do:

`nextflow run ikmb/exome-seq -r 1.5 --samples Samples.csv --assembly GRCh38 --kit xGen --email 'hello@gmail.com`

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
The following human genome assembly versions are supported on MedCluster (see resources.config on how this is set):

- GRCh37 (the 1000 genomes reference with decoys)
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

### `--baits` | `--targets` 
If you have used any other type of kit for your enrichment, you are able to provide the target and bait definitions from the command line during execution using `--baits` and 
`--targets`, respectively. Please note that these files must be in the Picard 
[interval_list](https://gatkforums.broadinstitute.org/gatk/discussion/1319/collected-faqs-about-interval-lists) format and have to be matched 
to the genome assembly (i.e. must have identical dictionary headers). 

### `--panel`
For practical reasons, it can be desirable to determine the coverage of a discrete set of target genes, such as for a gene panel. The pipeline currently 
supports the following panels:

- Dilatative Kardiomyopathie [cardio_dilatative]
- Hypertrophe Kardiomyopathie [cardio_hypertrophic]
- Non-Compaction Kardiomyopathie [cardio_non_compaction]
- Gene Immundefekt Agammaglobulinämie (25kb panel) [IMM_AGG]
- Gene Immundefekt Hypogammaglobulinämie (25kb panel) [IMM_HGG]
- Gene Immundefekt großes Panel [IMM]
- Gene Immundefekt intestinal (25kb panel) [IMM_IBD]
- Breast cancer panel [ breast_cancer ]
- Liver disease [ Liver ]
- Intellectual disability [ Intellectual_disability ]

Please note that this will also create additional run metrics, including a per-sample list of target exons that fall below a minimum sequence coverage. 

### `--all_panels`
This is a short-cut function to enable the production of statistics for all currently defined panels (for a given reference assembly!). Mutually exclusive with `--panel` and `--panel_intervals`. 

### `--panel_intervals`
This option allows the user to run non-defined panels. Must be in picard interval list format and match the sequence dictionary of the
genome assembly to run against (use with care!!!). Usually, you would start with a target list in BED format and convert this into an interval list
using the Picard Tools "BedToIntervalList" command.

### `--cnv`
Enable CNV calling using CNVkit. This option requires a pre-configured CNVkit reference matching the kit and assembly used for capture and mapping, respectively. Currently, this is only available for GRCh38 and xGen_v2. Alternatively, an external reference can be provided using the developer option `--cnv_ref`.

### `--vep`
Run variant effect prediction on the final VCF file(s). This option requires a locally available EnsEMBL cache and some databases (see cluster profiles for examples). 

## `--joint_calling` [ true (default) | false ]
Run joint calling on the samples rather than simply merging them down into one final VCF without generating sample-overarching genotyping for all possible sites. 

### `--kill`
For panel-based statistics, it is desirable to mark any exons that are known to underperform in exome sequencing - for example due to homology and
resulting multimapping (MAPQ = 0). This options allows the user to provide a list of panel targets that are to be listed as "KNOWN BAD" when compiling the
coverage report. An example is included for the IDT xGen v2 kit and assembly GRCh38 [here](../assets/kits/hg38_no_alt/idt_xgen_v2/kill.txt) .

### `--max_length` 
Set this to a positive number to trim all reads down to a desired size. Default: no size-trimming.

### `--interval_padding`
Set this to a positive number to include flanking regions of exon targets in the analysis. Default: 10

### `--skip_multiqc`
Skip the sending of a QC report. Default: false

### `--cram`
Create CRAM instead of BAM files to save space. Note that CRAM files are slower to read by IGV. 

### `--run_name`
Give this run a meaningful name (like a LIMS or project ID)

### `--fasta`
Provide path to a genome sequence in FASTA format (default: false, uses a pre-configured genome)

### `--dict`
Provide path to a genome sequence dictionary file (default: false, uses a pre-configured dictionary)

### `--dbsnp`
Provide path to a dbSNP reference VCF file for variant filtering (default: false, uses a pre-configured reference)

## Debug / custom arguments

### `--panel_coverage`
This option changes the cut-off for reporting lowly covered panel intervals (default: 10)

### `--cnv_ref`
Option to pass a custom CNVkit reference (cnn.gz) to the pipeline. This reference must match both the exom kit and assembly used! File must be compressed with gzip. 
