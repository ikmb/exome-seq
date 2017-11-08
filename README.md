![](images/ikmb_bfx_logo.png)

# IKMB Diagonistc Exome Pipeline

# Please note

This pipeline is under development. Specifically, the user can choose one of two processing chains for variant calling using either Freebayes or GATK4. The Freebayes workflow should work as expected. However, GATK4 is still being worked on by the BROAD Institute and the recommended processing chain not yet finalized. 

Both Freebayes and GATK4 are open-source - i.e. as user, you do not have to worry about licensing fees.

# Overview

This pipeline performs exome analysis on a set of samples. The following steps are included:

1. Read trimming (Trimmomatic 0.36)

2. Read alignment (BWA 0.7.15, Samtools 1.5)

3. Marking of duplicates (GATK4b5 or Picard 2.9.2)

4a. GATK4 workflow

* Base quality recalibration (GATK4b5)

* Variant calling using HaplotypeCaller (GATK4b5) 

* Merged genotype calls using GenotypeGVCFs (GATK4b5) [Includes 17 IKMB control exomes)

* Recalibration of SNPs and Indels (GATK4b5)

* Merged gVCF with control exomes removed

4b. FreeBayes workflow

* Joint Variant calling (Freebayes 1.1.0)

* Hard filtering of resulting VCF file (VCFtools)

5. Effect prediction with VEP (EnsEMBL86) and Annovar (2017)

## Installing the pipeline

To install this pipeline, simply clone the repository to a location on RZCluster:

`git clone git@git.ikmb.uni-kiel.de:bfx-core/NF-diagnostics-exome.git`

To update the code, run git update inside of the local clone:

`git update`

The pipeline is set up to work on RZCluster using the IKMB module system. Please make sure that you have set up this environment before launching a pipeline ru$

## Valid assemblies

The pipeline currently (officially) supports two genome assemblies - the "full" hg19 assembly from the GATK resource bundle and a stripped-down version that removes
any ALT contigs (hg19_clinical). To choose an assembly, use:

`--assembly hg19`

`--assembly hg19_clinical`

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

`ruby /path/to/git/samplesheet_from_folder.rb -f /path/to/fastq_folder`

where /path/to/fastq_folder points to a directory from which to read the fastq files. This requires a Ruby version >= 2.0 (available as software module on RZCluster.

WARNING: The script uses the base library ID to group read files; this will usually also take care of merging libraries that were sequenced across multiple lanes. However, it is advisable to make sure that each pair of FastQ files has the appropriate library and sample ID to accurately represent the experimental setup. 

### Read group ID tags

The read group ID and platform unit can be derived from the fastq headers like so:

Header: `@J00124:23:HGJJMBBXX:3:1101:31477:1138 1:N:0:CGTACTAG+ATAGAGAG`

RG_ID : `HGJJMBBXX.3.ATAGAGAG`

PU : `HGJJMBBXX.3` (This may not be 100% accurate, but available information is vague...)

## Executing the pipeline

The pipeline can be run as follows:

`nextflow -c /path/to/git/repo/nextflow.config run /path/to/git/repo/main.nf --samples /path/to/sample_list.csv`

Should the pipeline crash, you can try and resume it (after the problem has been fixed) adding "-resume" to the execution. 

### Supported Tool chains

This pipeline can be run with one of two processing chains - GATK4 and Freebayes. 

To use one of the two processing chains, use the `tool` option:

`nextflow run bfx-core/NF-diagnostics-exome --samples /path/to/sample_list.csv -hub ikmb --tool gatk`

`nextflow run bfx-core/NF-diagnostics-exome --samples /path/to/sample_list.csv -hub ikmb --tool freebayes`

Note that Freebayes is the default and does not need to be specified. 

### Supported enrichment kits

The pipeline is built to support more than one Exome kit. These can be selected from the command line using the --kit option.

`--kit Nextera` (the nextera rapid exome kit, 2017)

`--kit xGen` (the IDT xGEN panel, version 1.0)

This option defaults to "Nextera", so make sure this is in fact the kit you have used!

### Email notification

The pipeline will send out an Email notification upon completion if a recipient is speciefied via the flag `--email`:

`--email 'your.username@ikmb.uni-kiel.de'` 

## Output

By default, the pipeline output will be stored in the "output" subfolder from where you ran the nextflow process. You can provide an alternative location 
by specifying "--outdir /some/other/folder" on the command line. 

The output will be spread across several folders - one for each sample, one for the combined variant calls and a generic folder for library statistics and such. 

A brief description of the pipeline outputs is available under http://git.ikmb.uni-kiel.de/bfx-core/NF-diagnostics-exome/tree/master/doc
