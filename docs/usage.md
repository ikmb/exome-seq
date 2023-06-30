# Usage information

## Basic execution

The following command will execute this pipeline; the options will be discussed in the following.

If you are at the IKMB, the following will work:

`nextflow run ikmb/exome-seq --samples Samples.csv --assembly GRCh38 --kit xGen_v2 --email 'hello@gmail.com' --tools 'strelka,deepvariant,manta'` 

If you need to run a specific "release" of the pipeline, you can do:

`nextflow run ikmb/exome-seq -r 4.3 --samples Samples.csv --assembly GRCh38 --kit xGen_v2 --email 'hello@gmail.com` --tools 'strelka,deepvariant,manta'

If you try to run the pipeline on another system, you will need to configure a profile (see the installation instructions):

`nextflow run /path/to/main.nf --samples Samples.csv --assembly GRCh38 --kit xGen_v2 --tools 'deepvariant,expansionhunter,manta' --email 'hello@gmail.com' -profile your_profile`

## Mandatory arguments

### `--samples`
Information to this pipeline is given in form of a CSV sample sheet. This is to allow relevant, and possibly quite important, metadata to be included (such
as meaningful sample names, sequencing center etc).

The structure looks as follows:

```
patient;sample;status;library;readgroup;platform_unit;center;date;R1;R2
```
where status refers to the tumor status (0 = normal, 1 = tumor). The pipeline will automatically determine if somatic calling an be performed with or without a matched normal based on the sample and patient IDs and the tools that were requested. 

For convenience, we have included a simple script (bin/samplesheet_from_folder.rb) which accepts the path to a folder fill of PE exome data and automatically
writes a basic sample sheet. Obviously, it cannot derive meaningful sample and patient IDs from such information; this you would have to edit manually, if
you so choose.

```bash
ruby bin/samplesheet_from_folder.rb -f /path/to/foler > Samplesheet.csv`
```

### `--assembly` 
The following human genome assembly versions are supported on MedCluster (see resources.config on how this is set):

* GRCh38 (patch 1, with decoys and masked PAR regions - see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) )
* GRCh38_g1k (patch1, as used by the 1000 genomes consortium and the Illumina Dragen system - see [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/))
* GRCh38_no_alt (patch 1, no ALT contigs, with decoys and masked PAR regions - see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) )
* GRCh38_p14 (patch 14 without further modifications, see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/) )
* GRCh38_no_alt_p14 (patch 14 without ALT contigs, see [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/) )
* hg38 (the BROAD version of GRCh38, part of the GATK bundle, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) )

The choice is with the user, although we recommend a non-ALT version because downstream variant callers and effect prediction tools are not able to deal with ALT-located variants in a meaningful way. 

Basically: If you want a version that is optimized for short read alignment, use GRCh38_no_alt. If you need your data to be compatible with results from the Illumina Dragen platform, use GRCh38_g1k. If you need to work with the latest patch level, use GRCh38_no_alt_p14. 

### `--genomes_base` 
The root directory of the pre-installed indices. See "--build_references" on how to create this folder structure. 

### `--build_references`
This option is to be used instead of the primary workflow and will build the necessary index files for the respective reference assemblies. Please see our installation [instructions](installation.md) for details.

### `--aligner` [default = "bwa2"]
The following alignment algorithms are supported

- [BWA](https://github.com/lh3/bwa) (bwa)
- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) (bwa2)
- [DRAGMAP](https://github.com/Illumina/DRAGMAP) (dragmap)

Again, the choice is with the user. Note that BWA and BWA2 deliver near-identical results, but BWA2 is roughly 2-times faster. If you want to produce alignments equivalent to the Illumina Dragen system (v4), you will have to use the assembly 'GRCh38_g1k' and dragmap as your aligner. 

### `--tools`
The pipeline offers various tools for the analysis of variant information. Specifically:

1. SNPS and INDELs
   - [Deepvariant](https://github.com/google/deepvariant) (deepvariant)
   - [Strelka](https://github.com/Illumina/strelka) (strelka)
   - [Haplotyecaller](https://github.com/broadinstitute/gatk) (haplotypecaller)
2. Somatic variant calling
   - [Mutect2](https://github.com/broadinstitute/gatk) (mutect2)
2. Structural variants
   - [Manta](https://github.com/Illumina/manta) (manta)
3. Repeat expansions
   - [Expansion Hunter](https://github.com/Illumina/ExpansionHunter) (expansionhunter)
4. Copy number variants
   - [CNVkit](https://cnvkit.readthedocs.io/en/stable/) (cnvkit)
5. Protein-level haplotypes
   - [Haplosaurus](https://www.ensembl.org/info/docs/tools/vep/haplo/index.html) (haplosaurus)
   - [CSQ](https://samtools.github.io/bcftools/howtos/csq-calling.html) (csq)
6. HLA calling
   - [xHLA](https://github.com/humanlongevity/HLA) (xhla)
7. Variant effects:
   - [VEP](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)  (vep)

Your tools of choice can be provided like so:

`nextflow run ikmb/exome-seq --samples Samples.csv --assembly GRCh38 --tools 'deepvariant,manta,expansionhunter,strelka,cnvkit'`

If no tools are selected, the pipeline will stop after the deduplication of read alignments. 

Please also note that certain parts of the pipeline will be silently activated if a downstream tool requires them. Specifically, if you request variant calling with Strelka, the Manta pipeline for detecting structural variants will run whether you request it or not. This is because Strelka has been shown to perform better when informed by Manta INDEL predictions. 

### `--joint_calling`
The pipeline produces multi-vcf files through merging of the single-sample callsets. However, you can alternatively request for the samples to be called jointly, i.e. all in one process. This will cause individual callsets to take into consideration information from other samples, and result in numerous ref calls in the individual VCF files. An advantage of this approach is that individual sites can obtain support from multiple samples. This is typically done when analysing cohorts. Depending on the caller, the exact process of joint calling may differ. 

### `--kit`
Each exome capture kit has a target and a bait definition, i.e. information about the exons it enriches and the specific RNA bait sequences that are used 
for capture. This information is important so the pipeline knows which regions of the genome to analyze and how to compute run metrics. 

We have included these files for the following kits and genome assemblies:

`xGen_v2` (v2 release of the IDT xGen kit) [all assemblies]
`Agilent_v7` (v7 release of the Agilent SureSelect kit) [all assemblies]

### `--email`
Your Email address in quotes to which the pipeline report is sent upon completion. 

## Tool-specific options

### CNVkit
CNVKit runs in two separate modes: single and paired. In single mode, every patient sample will run stand-alone using a reference panel. This reference panel can be provided by the user through the option `--cnv_gz` (also see [here](https://cnvkit.readthedocs.io/en/stable/pipeline.html)) or will be built on-the-fly as a flat reference. Using a flat reference as usually not a good idea for exome data, since bait capture efficiency is not uniform, and neither is the coverage across exons. 
Paired mode on the other hand compares a normal sample with all tumor samples from the same patient and does not need a pre-existing reference panel.
Naturally, single mode always runs whereas the paired mode only activates if the right data is provided (i.e. at least one tumor and one normal sample from the same patient). 

For xGen v2, a reference panel is included and used automatically. Since this panel was produced from our own in-house data, you may want to override this option tho. 

#### `--cnv_gz`
If you wish to provide a CNVKit reference panel, this option is for you. This file must be compressed with gzip (.cnn.gz) and match the assembly and exome kit!
#### `--cnvkit_mode` [ default = "hmm-germline" ]
The segmentation mode for CNV intervals. Default is hmm-germline. Other options are documented [here](https://cnvkit.readthedocs.io/en/stable/pipeline.html#segment).
#### `--cnvkit_mode_tumor` [ default = "cbs" ]
The segmentation mode for CNV intervals when analysing tumor samples. Default is cbs. Other options are documented [here](https://cnvkit.readthedocs.io/en/stable/pipeline.html#segment).
### Deepvariant
#### `--glnexus_config` [ default = "DeepVariant" ]
The filter profile for gVCF merging in GLXNexus (DeepVariant). The default (DeepVariant) is fairly unconstrained. Other options are DeepVariantWGS and DeepVariantWES.
### GATK
#### `--gatk_hard_filter` [default = "ExcessHet > 54.69"]
This option allows users to specify on which annotations to hard-filter the GATK callset.
### MUTECT2
#### `--mutect_normals` [ default = null ]
Provide a matching panel of normals to help variant mutect2 variant filtration. Must be generated against the same reference assembly using the same capture kit (and preferably sequencing platform). For further instructions, see [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2).
#### `skip_mutect_pon` [ default = false ]
If your profile uses a pre-configued panel of normals, set this to true to skip it anyway.
### VEP
#### `--dbnsfp_db` [default = null]
Path to a copy of the dbNSFP database. 
#### `--dbscsnv_db` [default = null]
Path to a copy of the dbSCSNV database.
#### `--cadd_snps` [default = null]
Path to a local copy of the CADD SNPs in VCF format. 
#### `--cadd_indels` [default = null]
Path to a local copy of the CADD Indels in VCF format
#### `--vep_mastermind` default = null]
Path to a local copy of the Mastermind database.
## Calling regions and gene panels
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
### `--panel_coverage`
This option changes the cut-off for reporting lowly covered panel intervals (default: 10)
## Misc arguments

### `--kill`
For panel-based statistics, it is desirable to mark any exons that are known to underperform in exome sequencing - for example due to homology and
resulting multimapping (MAPQ = 0). This options allows the user to provide a list of panel targets that are to be listed as "KNOWN BAD" when compiling the
coverage report. An example is included for the IDT xGen v2 kit and assembly GRCh38 [here](../assets/kits/hg38_no_alt/idt_xgen_v2/kill.txt) .

### `--interval_padding`
Set this to a positive number to include flanking regions of exon targets in the analysis. Default: 10

### `--skip_multiqc`
Skip the sending of a QC report. Default: false

### `--run_name`
Give this run a meaningful name (like a LIMS or project ID)

## Expert options

### `--amplicon_bed`
A BED file specifying the location of amplicon primer positions. These will be masked from the final BAM file; no deduplication will be performed. Must match the assembly version. 

See the samtools documentation for an example of how this file needs to be formatted [here](http://www.htslib.org/doc/samtools-ampliconstats.html).

| Chromosome | Start | Stop | primer_name | score | strand |
| ---------- | ----- | ---- | ----------- | ----- | ------ |
| MN908947.3 | 1875  | 1897 | nCoV-2019_7_LEFT | 60 | + |
| MN908947.3 | 1868  | 1890 |  nCoV-2019_7_LEFT_alt0 | 60 | + |
| MN908947.3 | 2247  | 2269 | nCoV-2019_7_RIGHT | 60 | - |

