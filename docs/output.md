# Output

## Folder structure

All outputs, by default, will be written into the folder "results". Within this folder, you will find the following:

* One of more "patient" folders, named after the "IndividualID" provided through the samplesheet. Within each patient folder, there will be one folder per sample 
associated with that patient and therein the BAM files (and some statistical outputs).

* A folder "pipeline_info" containing basic information about this pipeline run, such as run time, the workflow graph etc. 

* A folder "Summary" containing the QC reports for various levels of data (fastq files, libraries, samples and, optionally, panels)

* A folder "Variants" containing the variant calls in VCF format, sorted by filtering procedure (HardFilter or VQSR)

## QC reports

### Fastq MultiQC

This QC report shows basic information on the individual fastq files - how much of the reads were classified as PCR duplicates, the mean GC content of the reads, percent of reads passing the length and
quality filters and finally how much of the read data can be assigned to sequencing adapters. 

For a good sample, you want relatively constant GC values across fastQ files, a high fraction of reads passing the filters (> 98%), and only little adapter content (< 5%). Irregular values
in any of these metrics may indicate issues with the sequencing or, more likely, the input material.

### Library MultiQC

The Library QC report looks at metrics relevant for the original sequencing library - specifically level of duplication. Duplication is a common issue in library protocols using some form of
selection/amplification. The level of duplication should be low, but can sometimes be fairly high - which means that you lose data for the downstream interpretation. 

### Sample MultiQC

The sample report contains essential information about the actual enrichment performance, including values critical for the final reporting, such as mean coverage of the exome target. 
This report furthermore contains a simple sex check by reporting the coverage of the SRY gene. Male samples sould have very high coverage (>10.000) in this region, 
wheras in female samples coverage should be very low (<<1000). If the value for any given sample is 0, it means that SRY was not part of the enrichment kit. 

### Panel MultiQC

This QC report contains a subset of the sample report, focused on a list of genes of interest for specific diagnostic applications. 



