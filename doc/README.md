![](images/ikmb_bfx_logo.png)

# IKMB Diagnostic Exome Pipeline -  Outputs

This pipeline performs the GATK best practice pipeline for variant calling (GATK4). Please not the various optional arguments to tweak the behavior of this pipeline. 
By default, this pipeline produces read alignments in BAM format and performs hard-filtering on the final variant list only. In order to also run variant score recalibration,
please use the --vsqr true option from the command line and note that this will only work with 30 or more exomes at once. 

## Reference data

If not indicated otherwise, analyses are run against the hg19 reference dataset distributed by the BROAD Institute. To obtain your own copy, please see:

https://software.broadinstitute.org/gatk/download/bundle

## Reporting

Located under 'output/Summary'

Reports are generated on the basis of the original FastQ files, the sequencing libraries and the (biological) samples. For each category, there
is one report summary in HTML format - use any up-to-date browser to visualize the results. 

The outputs should be fairly self-explanatory; note that the HTML outputs contain tool-tips in various places to help with the interpretation. 

## Read alignments (BAM or CRAM)

The primary output of this pipeline are the aligned, duplicate-marked reads. These can be found in the folder:

output/IndividualID/SamplesID/Processing/MarkDuplicates

Read alignments are stored in either BAM or lossless CRAM format. 

## Primary variant calls

The primary VCF file will be located under:

output/Variants/Final/HardFiltered/some_file.vcf

## Annotated variant calls

VCFs may optionally be annotated using the EnsEMBL variant effect predictor (https://www.ensembl.org/info/docs/tools/vep/index.html). The annotated output will be located under output/Annotation/VEP/




