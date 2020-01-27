# Output

## Folder structure

All outputs, by default, will be written into the folder "results". Within this folder, you will find the following:

* One of more "patient" folders, named after the "IndividualID" provided through the samplesheet. Within each patient folder, there will be one folder per sample associated with that patient and therein the BAM files (and some statistical outputs).

* A folder "pipeline_info" containing basic information about this pipeline run, such as run time, the workflow graph etc. 

* A folder "Summary" containing the QC reports for various levels of data (fastq files, libraries, samples and, optionally, panels)

* A folder "Variants" containing the variant calls in VCF format, sorted by filtering procedure (HardFilter or VQSR)


