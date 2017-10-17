![](images/ikmb_bfx_logo.png)

# IKMB Diagnostic Exome Pipeline -  Outputs

This pipeline runs with one of two processing chains using either GATK4 or Freebayes 1.1.0. The output is structured as follows:

## Reference data

If not indicated otherwise, analyses are run against the hg19 reference dataset distributed by the BROAD Institute - minus ALT contigs (to avoid issues arising from mapping ambiguities in the presence of alternative haplotypes). To obtain your own copy, please see:

https://software.broadinstitute.org/gatk/download/bundle

## Reporting

Located under 'output/Summary'

Reports are generated on the basis of the original FastQ files, the sequencing libraries and the (biological) samples. For each category, there
is one report summary in HTML format - use any up-to-date browser to visualize the results. 

The outputs should be fairly self-explanatory; note that the HTML outputs contain tool-tips in various places to help with the interpretation. 

## Finished BAM files

The primary output of this pipeline are the aligned, duplicate-marked reads. These can be found in the folder:

output/IndividualID/SamplesID/Processing/MarkDuplicates

## Primary variant calls

Variants are called using either Freebayes or GATK. The primary VCF file will be located under:

output/Variants/VariantCaller/some_file.vcf

## Annotated variant calls

The final VCF file (after filtering) is usually annotated using both Annovar (Mid 2017) and the EnsEMBL VEP (release 90). The resulting annotated VCF files can be found under

output/Annotation/AnnotationTool/some_file.vcf





