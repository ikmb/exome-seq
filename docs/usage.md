# Usage information

## Basic execution

The following command will execute this pipeline; the options will be discussed in the following:

`nextflow run main.nf --samples Samples.csv --assembly hg19 --kit xGen --email 'hello@gmail.com'`

Let's dissect this in the following:

### The samplesheet

Information to this pipeline is given in form of a CSV sample sheet. This is to allow relevant, and possibly quite important, metadata to be included (such as meaningful sample names, sequencing center etc). 

The basic format of the sample sheet is as follows:

`IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2`

That is to say:

- ID for the patient/probant that was sequenced
- ID of the sample from that patient/probant (such as a tumor or a "normal" sample)
- ID of the sequencing library
- read group ID (a unique ID for this piece of data, usually a combination of flow cell id and library id)
- platform unit (similar to the above)
- platform model - the type of sqeuencer this data was generated on (e.g. Novaseq 6000)
- Center - Name of the sequencing center
- Date - When this data was generated
- R1 - path to the left file of a paired-end data set
- R2 - patht to the right file of a paried-end data set

For convenience, we have included a simple script (bin/samplesheet_from_folder.rb) which accepts the path to a folder fill of PE exome data and automatically write a basic sample sheet. Obviously, it cannot derive meaningful sample and patient IDs from such information; this you would have to edit manually, if you so choose. 

`ruby bin/samplesheet_from_folder.rb -f /path/to/foler > Samplesheet.csv`

## The genome assembly

The choice of genome assembly depends a bit on your local community and mostly personal preference. If you downloaded the GATk genome bundle as per the installation instruations, you have the following options automatically available:

- GRCh37 (the 1000 genomes reference)
- GRCh38 (the current human reference assembly)
- hg19 (another version of GRCh37, also referred to as the UCSC reference)

## The exome kit

Each exome capture kit has a target and a bait definition, i.e. information about the exons it enriches and the specific RNA bait sequences that are used for capture. This information is important so the pipeline knows whichs regions of the genome to analyze and how to compute run metrics. 

We have included these files for two capture kits - IDT xGen and Nextera. YOu can choose one or the other:

`--kit Nextera`

`--kit xGen`

We also offer a custom version of xGen (`--kit xGen_custom`), which includes a few additional SNPs missing from xGen; this is really only relevant if you use our custom mix - so it's safe to ignore. 

## Reporting

The pipeline will send a basic report upon completion to your Email address of choice. This requires for the compute running the nextflow process to have a configured Mail demon.






