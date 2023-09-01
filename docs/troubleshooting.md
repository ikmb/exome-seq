# Troubleshooting

Any common issues we encounter will be listed here with possible solutions. 

## IGV is slow when loading CRAM files

This is a known limitation of IGV and can easily be solved by not using the built-in reference but instead load the fasta (+.fai) file via "Genomes>Load genome from file". This approach is recommended also because the IGV reference for hg38/GRCh38 is likely not the same as the one that was used for mapping and will hence throw validation errors in regions were the md5 sum of the mapping reference and built-in IGV reference differ. 

If you do so, you either have to load a matching annotation from the IGV data hub (File>Load from server) or provide it as GFF3 formatted file. 

## The pipeline crashes when I request t2t-2.0 as assembly

Support for t2t-2.0 is currently very experimental and some parts of the pipeline are missing the necessary reference files (because they have not yet been made available). 

Specifically, the following limitations apply:

* Aligner: Dragmap is not recommended for use with t2t because we have not yet investigated the proper build strategy for this combination (specifically re: masked regions etc). 
* Callers: GATK cannot be used, period. GATK requires several callibration references that are not available for t2t - and lift-over to convert the existing ones seems like a bad idea for this particular use case (i.e. calibration). 
* Exome kits: xGen_v2 is the only supported exome kit and is missing over 1000 targets due to lift-over errors (mostly on chrY)



