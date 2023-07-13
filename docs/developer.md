# Development guide

## Code structure

This pipelines uses Nextflow [DSL2](https://www.nextflow.io/docs/latest/dsl2.html) code conventions. The general structure can be broken down as follows:

- `main.nf` is the entry into the pipeline. It loads some library files and calls the actual workflow in `workflows/dragen.nf`
- `workflows/exome-seq.nf` defines the main logic of the pipeline. It sets some key options, reads the samplesheet and calls the various subworkflows and modules
- `subworkflows/` location where self-contained processing chains are defined (also part of the pipeline logic).
- `modules/` the various process definitions that make up the pipeline
- `conf/` location of the general and site-specific config files
- `conf/resources.config` holds many of the important options about the location of reference files etc
- `assets` holds some of the (smaller) reference files needed by the pipeline, such as exome kit interval lists
- `bin` location of custom scripts needed by some of the pipeline processes
- `doc` the documentation lives here

## Pipeline logic

For a very broad overview of what this pipeline does, you can check out the simplified [pipeline graph](pipeline.md). A more fine grained description follows below. 

### Summary

The ikmb/exome-seq pipeline performs variant detection on Illumina short read data that was captured using a specific exome capture kit. These kits specifically enrich coding regions of the 
genome to minimize the amount of sequencing information that is required to cover the most relevant parts (i.e. coding DNA) of the genome for genetic analysis. 

The basic workflow implements the following:

- Read trimming and quality control
- Alignment against a chosen reference assembly version
- Cleanup (sorting, deduplication, repairing invalid metadata) of the alignment
- Variant calling with one or several callers
  - Single nucleotide polymorphisms and insertions/deletions (SNPs, Indels)
  - Copy number variants
  - Structural variants
- Effect prediction of variant calls (mostly SNPs/Indels)
- Summary report(s)

This basic outline is defined in the primary workflow [exome-seq.nf](../workflows/exome-seq.nf). The actual processes are themselves included in the various [subworkflows](../subworkflows/). 

### Important details

#### Samplesheet

An important initial step is the validation of the input samplesheet through [VALIDATE_SAMPLESHEET](../modules/validate_samplesheet.nf). Any changes to the samplesheet have to 
added here to make sure that no error is thrown. Also make sure to check [below](#read-and-metadata-channel).

#### Read and metadata channel

While the input to the pipeline is a csv-formatted samplesheet, most channels pass a hash-map `meta` as well as any of the file objects that are emitted upstream. This hash-map is created inside the [TRIM_AND_ALIGN](../subworkflows/align.nf) subworkflow. If changes are made to the samplesheet, updates to this subworkflow are needed!

#### Normal versus somatic samples

This pipeline deals differently with input data depending on whether or not it comes from a normal (germline) sample and/or somatic sample. And even more complexity is added if a a patient as a normal and multiple somatic samples. The logic behind this comes into play both for GATK-specific processing as well as non-GATK processing. This is because GATK requires score-recalibrated BAM files (because: reasons) and all other tools do not. In either case, the logic is the same - tag BAM files into "normal" and "tumor", group Channels by patient and sample IDs and then check which of these groupings has a) no somatic samples, b) no normal samples or c) normal and one or more somatic samples. These resulting channels are then plugged into the appropriate processing chains. 




