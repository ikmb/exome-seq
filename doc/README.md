![](images/ikmb_bfx_logo.png)

# IKMB Diagnostic Exome Pipeline -  Outputs

This pipeline performs the GATK best practice pipeline for variant calling (GATK4). 

## Reference data

If not indicated otherwise, analyses are run against the hg19 reference dataset distributed by the BROAD Institute. To obtain your own copy, please see:

https://software.broadinstitute.org/gatk/download/bundle

## Reporting

Located under 'output/Summary'

Reports are generated on the basis of the original FastQ files, the sequencing libraries and the (biological) samples. For each category, there
is one report summary in HTML format - use any up-to-date browser to visualize the results. 

The outputs should be fairly self-explanatory; note that the HTML outputs contain tool-tips in various places to help with the interpretation. 

## Finished CRAM files

The primary output of this pipeline are the aligned, duplicate-marked reads. These can be found in the folder:

output/IndividualID/SamplesID/Processing/MarkDuplicates

Read alignments are stored in lossless CRAM format. 

## Primary variant calls

The primary VCF file will be located under:

output/Variants/Final/some_file.vcf

## Annotated variant calls

VCFs may optionally be annotated using the EnsEMBL variant effect predictor (https://www.ensembl.org/info/docs/tools/vep/index.html). The annotated output will be located under output/Annotation/VEP/




