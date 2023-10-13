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

This pipeline deals differently with input data depending on whether or not it comes from a normal (germline) sample and/or somatic sample. And even more complexity is added if a a patient as a normal and multiple somatic samples. 
The logic behind this comes into play both for GATK-specific processing as well as non-GATK processing. This is because GATK requires score-recalibrated BAM files (because: reasons) and all other tools do not. In either case, the logic is the same - tag BAM files into "normal" and "tumor", group Channels by patient and sample IDs and then check which of these groupings has a) no somatic samples, b) no normal samples or c) normal and one or more somatic samples. These 
resulting channels are then plugged into the appropriate processing chains. 

### Adding a new exome kit

Exome kits are used to limit the regions of interest during variant calling and are needed for various quality control metrics. Exome kits are added to this code base for [each assembly](https://github.com/ikmb/exome-seq/tree/master/assets/kits) separately. 

An exome kit consists of two files - one for the baits and one for the actual targets. Baits in this context are the locations of the actual RNA baits used during capture, whereas targets are those regions which the kit aims to make accessible for
analysis. Most vendors ship this information as two separate files. If however no bait file is provided, you can use the target file as bait file (i.e. copy it). However, this is not ideal and will likely mess with some of the metrics that
require bait intervals. 

The pipeline expects targets and baits to be present in [interval list](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) format. Interval lists are similar to BED files, but include a sequence
dictionary. This was introduced by GATK to enable validation of calling regions against assemblies to avoid accidential mismatches (presumably). For all non-GATK parts of the pipeline, interval lists will be converted to BED files on the fly. 

Exome kits are configured in the [resources.config](https://github.com/ikmb/exome-seq/blob/master/conf/resources.config) file. Here, each assembly has a hash key with all the associated files. You can have a look at one of the existing kits to figure
out how this works - and then add any new kit accordingly. It is recommended to add new kits to all the assemblies equally to avoid confusion over which kit is available for which assembly, so you'll have to write a few bash loops to 
[convert](https://gatk.broadinstitute.org/hc/en-us/articles/360037056472-BedToIntervalList-Picard-) your original BED file(s) to interval lists for each assembly. 

### Adding a new panel

Panels are specific lists of genes that are used to evaluate the success of exome-sequencing with respect to a specific kind of diagnostic application (i.e., different diagnostic approaches focus on different genes - and you may wish to 
know if your genes of interest were all adequately covered by the analysis). This is a left-over of when exome-seq was used in routine diagnostics and probably is of no further interest. If you, occasionally, need to evaluate the performance
of your analysis against a specific set of genes, consider using `--panel_intervals` to pass a gene panel in [interval list](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) format rather than adding it permanetly to this code base. 

If you are set on adding a new panel regardless, then here is a short how-to:

First, you need to obtain your genic regions of interest in BED format. Typically, these will contain the position of (coding) exons - one entry per exon. There exists a small EnsEMLB API [script](https://github.com/ikmb/exome-seq/blob/master/bin/ensembl_panel2bed.pl) 
that will produce this file for you, given a text file with canonical (HGNC-compliant) gene names. For this script to work, you must have a working installation of the [EnsEMBL Perl API](https://www.ensembl.org/info/docs/api/api_installation.html) 
installed. Note that the Perl API is version-locked. So to get data from a specific release of EnsEMBL, you must install that version specifically. 

Second, once you have your BED file, you should convert it into an [interval list]((https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists) and store it in the respective asset [folders](https://github.com/ikmb/exome-seq/tree/master/assets/panels).
Analogous to the exome kits, it is recommended to just set up a gene panel for all possible genome references to avoid unnecessary confusion. 

Finally, the pipeline needs to know about the new panel, which you can configure in [resources.config](https://github.com/ikmb/exome-seq/blob/master/conf/resources.config). Take any of the existing panels to see how this works. 

### Building a CNVkit reference

CNV calling is difficult to optimize for exome analysis, because exome data is subject to various biases - including capture efficency of the individual target regions. The resulting fluctuation in coverage interferes with the detection power of CNVkit, which
is based primarily on detecting divergence in coverage from an expected "normal". To address this issue, CNVkit accepts a panel of reference samples to learn what the normal distribution of coverage over individual targets looks like. 

However, such a reference set must meet certain criteria. Specifically, it must match the sample to be analyzed exactly, meaning it must have been produced with exactly the same processing chain as your sample of interest. This includes not only
the laboratory workflow, but also the exome kit and the sequencing platform. Any change in one of these parts will degrade the detection power of your CNV analysis.

So, what do you need?

It is recommended to collect at least 25-50 samples. These should all meet highest quality standards (coverage, read quality, etc) and must not be biased towards any given disease/phenotype to avoid training real CNVs as base line. 

Once this is in place, there is a detailed instruction on how to build a reference available from the [CNVkit documentation](https://cnvkit.readthedocs.io/en/stable/pipeline.html).

### Building a Mutect2 reference

Mutect2 requires a reference panel of normals to learn about baseline noise in your sequencing data. Similar to CNVkit, this reference panel must be matched very specifically to the data that you plan to analyze (lab workflow, exome kit and sequencer). 

To build such a panel, please follow the instructions provided by the developers [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132). 

As input, use a reasonably sized set of reference samples from "healthy" controls - at least 30-40, if possible. Make sure you produce the data with this pipeline, with the exome kit of interest and against your preferred reference assembly. 
Since you will be working with BAM files during the process of panel building, also make sure to enable gatk in the list of tools to run so that you get a base-score-recalibrated BAM for this purpose (by default, the pipeline will not perform 
this recalibration as it is only needed by GATK). 

``` 
nextflow run ikmb/exome-seq --samples samples.csv --assembly GRCh38_p14_no_alt --tools 'gatk' --kit xGen_v2
```

Once your have your BAM files, process each with Mutect2 to produce a VCF file:

```
 gatk Mutect2 -R /path/to/genome.fa -I sample.bam  --max-mnp-distance 0 -O sample.vcf.gz -OVI 
```

Next, import all your samples into a GATK genome db

```
gatk GenomicsDBImport -L exome_kit.interval_list --genomicsdb-workspace-path GenomeDB -V sample1.vcf.gz -V sample2.vcf.gz -V samplen.vcf.gz
```

Finally, use this genome db to build the panel of normals

```
CreateSomaticPanelOfNormals \
   -R /path/to/genome.fa \
   --germline-resource /work_ifs/ikmb_repository/references/gatk/v2/hg38_no_alt/af-only-gnomad.hg38.vcf.gz \
   -V gendb://GenomeDB \
   -O my_pon.vcf.gz \
   -OVI
```

Note that the germline resource used here is provided by GATK (see instructions above) and consists essentally of a stripped-down version of the gnomad variant reference. This file is already
present on the MedCluster, so can be grabbed from the location as shown in the example. Else, Google is your friend. 

Your Mutect2 panel of normals can then be added as a permanent resource to your site-specific [config file](https://github.com/ikmb/exome-seq/blob/master/conf/medcluster.config) - or you can pass it to the pipeline dynamically as
`--mutect_normals`. If you use the latter option, just be careful not to accidentially mismatch exome kit, assembly and panel of normals.  
	        	        	        	        	        

	                        	                        	                        









