# Troubleshooting

Any common issues we encounter will be listed here with possible solutions. 

## IGV is slow when loading CRAM files

This is a known limitation of IGV and can easily be solved by not using the built-in reference but instead load the fasta (+.fai) file via "Genomes>Load genome from file". This approach is recommended also because the IGV reference for hg38/GRCh38 is likely not the same as the one that was used for mapping and will hence throw validation errors in regions were the md5 sum of the mapping reference and built-in IGV reference differ. 

If you do so, you either have to load a matching annotation from the IGV data hub (File>Load from server) or provide it as GFF3 formatted file. 

