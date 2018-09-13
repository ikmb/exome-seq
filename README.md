![](images/ikmb_bfx_logo.png)

# IKMB Diagnostic Exome Pipeline

# Please note

This pipeline offers a end-to-end workflow for exome analysis using the GATK4 toolchain

- trimming with TrimGalore
- read alignment with BWA
- duplicate marking using Picard MarkDuplicates
- quality score recalibration
- gvcf calling
- joint variant calling
- variant recalibration (SNPs and Indels) and filtering
- effect prediction using VEP (optional)

## Installing the pipeline

To install this pipeline, simply clone the repository to a location on RZCluster:

`git clone git@git.ikmb.uni-kiel.de:bfx-core/NF-diagnostics-exome.git`

To update the code, run git update inside of the local clone:

`git update`

The pipeline is set up to work on RZCluster using the IKMB module system. Please make sure that you have set up this environment before launching a pipeline ru$

### Prerequisites

The following steps/resources are needed to run this pipelines:
* A working and configured Conda/Bioconda installation (please also see: https://bioconda.github.io/ - specifically point 2.) Set up channels)
* The GATK resource bundle (or equivalent resources). For a way to set up you own resources, please have a look at config/rzcluster.config)
* Bait/Target files matching your genome assembly of choice (files for assembly hg19 are included with this code base and can be used via the `--kit` flag)

## Valid assemblies

The pipeline currently (officially) supports the following genome assemblies from the GATK bundle:

`--assembly hg19` 
`--assembly hg38`
`--assembly b37`

Please note that only the exome kit files for hg19 are currently available; to use the other assemblies, you must specify your custom bait files from the command line using:

`--targets path/to/targets.interval_list`
`--baits path/to/baits.interval_list`

These files must be in the Picard "interval list" format and include matching sequence dictionaries to the genome sequence you wish to use them with. 

## Input format

The pipeline parses a config file to find the input files and relevant meta data. The expected format is semicolon-separated CSV with the following headers:

  * INDIVIDUAl_ID - The ID of the individual from which the sample was derived.
  * SAMPLE_ID - The ID of the sample. Note that more than one sample can come from the same individual (e.g. tumor/normal pair)
  * LIBRARY_ID - The ID of the DNA library. Multiple sequencing libraries can be prepared from the same sample.
  * RG_ID - Read group ID
  * PLATFORM_UNIT - Generally, this is the read group ID plus the library ID
  * PLATFORM - Sequencer (e.g. illumina)
  * PLATFORM_MODEL - Sequencer (e.g. HiSeq2500)
  * RUN_DATE - Date of sequencing run
  * CENTER - Location of sequencing run
  * R1 - Full path to Fastq file 1
  * R2 - Full path to Fastq file 2

A template file (SampleTemplate.csv) is included with this code base. 

You can use the script bin/samplesheet_from_folder.rb to generate a automatic sample sheet:

`ruby /path/to/git/bin/samplesheet_from_folder.rb -f /path/to/fastq_folder`

where /path/to/fastq_folder points to a directory from which to read the fastq files. This requires a Ruby version >= 2.0 (available as software module on RZCluster.

WARNING: The script uses the base library ID to group read files; this will usually also take care of merging libraries that were sequenced across multiple lanes. However, it is advisable to make sure that each pair of FastQ files has the appropriate library and sample ID to accurately represent the experimental setup. 

### Read group ID tags

The read group ID and platform unit can be derived from the fastq headers like so:

Header: `@J00124:23:HGJJMBBXX:3:1101:31477:1138 1:N:0:CGTACTAG+ATAGAGAG`

RG_ID : `HGJJMBBXX.3.ATAGAGAG`

PU : `HGJJMBBXX.3` (This may not be 100% accurate, but available information is vague...)

## Executing the pipeline

The pipeline can be run as follows:

This pipeline requires Java 1.8, Nextflow 0.31 or greater and conda/miniconda to run. On RZCluster, you can do:

`module load IKMB Java Nextflow miniconda2`

The command call for the pipeline itself:

`nextflow -c /path/to/git/repo/nextflow.config run /path/to/git/repo/main.nf --samples /path/to/sample_list.csv`

Should the pipeline crash, you can try and resume it (after the problem has been fixed) adding "-resume" to the execution. 

### Supported enrichment kits

The pipeline is built to support more than one Exome kit. These can be selected from the command line using the --kit option.

`--kit Nextera` (the nextera rapid exome kit, 2017)

`--kit xGen` (the IDT xGEN panel, version 1.0)

`--kit xGen_custom`(the IDT xGEN panel, version 1.0, with additional custom targets)

This option defaults to "Nextera", so make sure this is in fact the kit you have used!

### Email notification

The pipeline will send out an Email notification upon completion if a recipient is speciefied via the flag `--email`:

`--email 'your.username@ikmb.uni-kiel.de'` 

## Output

By default, the pipeline output will be stored in the "output" subfolder from where you ran the nextflow process. You can provide an alternative location 
by specifying "--outdir /some/other/folder" on the command line. 

Within the output folder will be three subfolders:

- Common (files common to all tool chains - i.e. duplicate marked read alignments and alignment statistics)
- Variants (the joint, filtered variant calls)
- Summary - graphical summary reports for Fastq files, libraries and samples
- Individual data (finalized read alignments, alignment statistics etc)
