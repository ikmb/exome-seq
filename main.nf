#!/usr/bin/env nextflow

/**
===============================
IKMB Diagnostic Exome Pipeline
===============================

This Pipeline performs one of two workflows to generate variant calls and effect predictions
using either the GATK processing chain or Freebayes.

### Homepage / git
http://git.ikmb.uni-kiel.de/bfx-core/NF-diagnostics-exome
### Implementation
Implemented in Q1 2018

This pipeline is based on the updated GATK best-practices (where applicable).
 - trimming (Trimgalore)
 - Alignment (BWA) and dedup
 - variant calling (Strelka)
 - variant effect prediction (VEP)

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version
VERSION = "1.0-alpha1"
params.version = VERSION

// Help message
helpMessage = """
===============================================================================
IKMB Diagnostic Exome pipeline | version ${VERSION}
===============================================================================
Usage: nextflow -c /path/to/git/nextflow.config run /path/to/git/main.nf --assembly hg19_clinical --kit Nextera --samples Samples.csv
This example will perform an exome analysis against the hg19 (with decoys) assembly, assuming that exome reads were generated with
the Nextera kit and using the GATK4 best-practice workflow. 
Required parameters:
--samples                      A sample list in CSV format (see website for formatting hints)
--assembly                     Name of the reference assembly to use
--effect_prediction	       Whether to run effect prediction on the final variant set (default: false)
Optional parameters:
--run_name 		       A descriptive name for this pipeline run
--fasta				A reference genome in FASTA format (set automatically if using --assembly)
--dbsnp				dbSNP data in VCF format (set automatically if using --assembly)
Output:
--outdir                       Local directory to which all output is written (default: output)
Exome kit:
--kit                          Exome kit used (default: Nextera)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// #############
// INPUT OPTIONS
// #############

// Sample input file
inputFile = file(params.samples)

// Giving this pipeline run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "hg19"

REF = params.fasta ?: file(params.genomes[ params.assembly ].fasta)
DBSNP = params.dbsnp ?: file(params.genomes[ params.assembly ].dbsnp )
VEP_CACHE = params.vep_cache

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits
TARGET_BED = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets_bed

params.effect_prediction = true
params.hard_filter = false

// Location of applications used
OUTDIR = file(params.outdir)

// Available exome kits

// Whether to send a notification upon workflow completion
params.email = false

// Whether to use a local scratch disc
use_scratch = params.scratch

// Make sure the Nextflow version is current enough
try {
    if( ! nextflow.version.matches(">= $params.nextflow_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nextflow_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please use a more recent version of Nextflow!\n" +
              "============================================================"
}

logParams(params, "pipeline_parameters.txt")

// Header log info
log.info "========================================="
log.info "IKMB Diagnostic Exome pipeline v${VERSION}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${params.run_name}"
log.info "========================================="

// Read sample file 
Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  readPairsFastp }

process runFastp {

	tag "${indivID}|${sampleID}"

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, date, fastqR1, fastqR2 from readPairsFastp

	output:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, date, center, file(left),file(right) into inputBwa
	set file(html),file(json) into fastp_results

	script:
	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right -w ${task.cpus} -j $json -h $html --length_required 35 --cut_by_quality3
	"""
}

process runBWA {

    tag "${indivID}|${sampleID}|${libraryID}|${rgID}"
    // publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/", mode: 'copy'

    //scratch use_scratch
	
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
    output:
    set indivID, sampleID, file(outfile) into runBWAOutput
    
    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"	

    """
	bwa mem -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $left $right | samtools sort -O bam -m 2G -@ 4 - > $outfile
    """	
}

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])

process mergeBamFiles_bySample {

        tag "${indivID}|${sampleID}"
	
	input:
        set indivID, sampleID, file(aligned_bam_list) from runBWAOutput_grouped_by_sample

	output:
	set indivID,sampleID,file(merged_bam) into mergedBamFile_by_Sample

	script:
	merged_bam = sampleID + ".merged.bam"

	"""
		picard MergeSamFiles \
			INPUT=${aligned_bam_list.join(' INPUT=')} \
			OUTPUT=merged.bam \
			CREATE_INDEX=false \
			CREATE_MD5_FILE=false \
			SORT_ORDER=coordinate

		picard SetNmMdAndUqTags \
			REFERENCE_SEQUENCE=$REF \
			INPUT=merged.bam \
			IS_BISULFITE_SEQUENCE=false \
			OUTPUT=${merged_bam}
	"""
}

process runMarkDuplicates {

	tag "${indivID}|${sampleID}"
        publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

        // scratch use_scratch

        input:
        set indivID, sampleID, file(merged_bam) from mergedBamFile_by_Sample

        output:
        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into MarkDuplicatesOutput, BamForMultipleMetrics, runPrintReadsOutput_for_OxoG_Metrics, runPrintReadsOutput_for_HC_Metrics, BamForDepthOfCoverage
	file(outfile_md5) into MarkDuplicatesMD5
	file(outfile_metrics) into DuplicatesOutput_QC
	file(outfile_bam) into (inputStrelka,inputManta)
	file(outfile_bai) into (inputStrelkaBai,inputMantaBai)

        script:
        outfile_bam = sampleID + ".dedup.bam"
        outfile_bai = sampleID + ".dedup.bai"
	outfile_md5 = sampleID + ".dedup.bam.md5"

        outfile_metrics = sampleID + "_duplicate_metrics.txt"

	"""
        	gatk --java-options "-Xms4G -Xmx${task.memory.toGiga()-1}G" MarkDuplicates \
                	-I ${merged_bam} \
	                -O ${outfile_bam} \
        	        -M ${outfile_metrics} \
                        --CREATE_INDEX true \
			--ASSUME_SORT_ORDER=coordinate \
			--MAX_RECORDS_IN_RAM 100000 \
			--CREATE_MD5_FILE true \
                        --TMP_DIR tmp
	"""

}

// ------------------------------------------------------------------------------------------------------------
//
// Perform base quality score recalibration (BQSR) including
// 1) Generate a recalibration table
// 2) Generate a new table after applying recalibration
// 3) Compare differences between recalibration tables
// 4) Apply recalibration
//
// ------------------------------------------------------------------------------------------------------------

process runManta {


	tag "ALL"
        publishDir "${OUTDIR}/Manta/Variants", mode: 'copy'

	input:
	file(bams) from inputManta.collect()
	file(indices) from inputMantaBai.collect()

	output:
	file("*.vcf.gz") into outputManta
	file("candidateSmallIndels.vcf.gz") into mantaCandidates

	script:
	"""

	configManta.py \
	-bam=${bams.join(' --bam=')} \
        --referenceFasta $REF \
        --exome \
        --callRegions $TARGET_BED \
        --runDir Manta

	python Manta/runWorkflow.py -m local -j ${task.cpus}

	cp Manta/results/variants/*.vcf.gz* .

	"""
}


process runStrelka {

        tag "ALL"
	publishDir "${OUTDIR}/Strelka/Variants", mode: 'copy'

	input:
	file(bams) from inputStrelka.collect()
	file(indices) from inputStrelkaBai.collect()

	output:
	set file(vcf), file(vcf_index) into outputStrelka

	script:

	vcf = "${params.run_name}.variants.vcf.gz"
	vcf_index = vcf + ".tbi"

	"""
		configureStrelkaGermlineWorkflow.py \
		--bam=${bams.join(' --bam=')} \
		--referenceFasta $REF \
		--exome \
		--callRegions $TARGET_BED \
		--runDir Strelka

		python Strelka/runWorkflow.py -m local -j ${task.cpus}

		bcftools annotate -a $DBSNP -c ID Strelka/results/variants/variants.vcf.gz | bgzip -c > $vcf
		tabix $vcf
		
	"""
}


process runSplitStrelkaVcf {
	
	tag "ALL"
	publishDir "${OUTDIR}/Strelka/Variants/BySample", mode: 'copy'

	input:
	set file(vcf),file(index) from outputStrelka
	
	output:
	set file("*.vcf.gz"),file("*.vcf.gz.tbi") into inputVep

	script:

	"""
		for sample in `bcftools query -l $vcf`; do bcftools view -f "PASS" -s \$sample $vcf | bcftools filter -i 'GT!="./."' -i 'GT!="."' -i 'GT!="0/0"' | python $baseDir/bin/filter_strelka_vcf.py | bgzip -c > \$sample.vcf.gz && tabix \$sample.vcf.gz ; done;
	"""

}

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process runCollectMultipleMetrics {
	tag "${indivID}|${sampleID}"
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'
 
	scratch use_scratch
	    
	input:
	set indivID, sampleID, bam, bai from BamForMultipleMetrics

	output:
	file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

	script:       
	prefix = sampleID + "."

	"""
		picard -Xmx5g CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
       	        PROGRAM=CollectSequencingArtifactMetrics \
                PROGRAM=CollectQualityYieldMetrics \
	        PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectBaseDistributionByCycle \
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${REF} \
		DB_SNP=${DBSNP} \
		INTERVALS=${BAITS} \
		ASSUME_SORTED=true \
		QUIET=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	

process runHybridCaptureMetrics {

    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_HC_Metrics

    output:
    file(outfile) into HybridCaptureMetricsOutput mode flatten

    script:
    outfile = sampleID + ".hybrid_selection_metrics.txt"

    """
        picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                TARGET_INTERVALS=${TARGETS} \
                BAIT_INTERVALS=${BAITS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

process runOxoGMetrics {

    tag "${indivID}|${sampleID}"
    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_OxoG_Metrics

    output:
    file(outfile) into runOxoGMetricsOutput mode flatten

    script:
    outfile = sampleID + ".OxoG_metrics.txt"

    """

         picard -Xmx${task.memory.toGiga()}G CollectOxoGMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                DB_SNP=${DBSNP} \
                INTERVALS=${TARGETS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

// ------------------------------------------------------------------------------------------------------------
//
// Plot results with multiqc
//
// ------------------------------------------------------------------------------------------------------------

process runMultiqcFastq {

    tag "Generating fastq level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Fastq", mode: 'copy'
	    
    input:
    file('*') from fastp_results.flatten().toList()
    
    output:
    file("fastp_multiqc*") into runMultiQCFastqOutput
    	
    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n fastp_multiqc *.json *.html
    """
}

process runMultiqcLibrary {

    tag "Generating library level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Library", mode: 'copy'
	    
    input:
    file('*') from DuplicatesOutput_QC.flatten().toList()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput
    	
    script:

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *.txt
    """
}

process runMultiqcSample {

    tag "Generating sample level summary and QC plots"
    publishDir "${OUTDIR}/Summary/Sample", mode: 'copy'
	    
    input:
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file('*') from runOxoGMetricsOutput.flatten().toList()
        
    output:
    file("sample_multiqc.html") into runMultiQCSampleOutput
    	
    script:

    def subject = 'Diagnostic exome analysis quality report'
    def recipient = params.email

    """
    cp $baseDir/config/multiqc_config.yaml multiqc_config.yaml
    multiqc -n sample_multiqc *

    """
}

// *************************
// Variant effect prediction
// *************************

process runVep {

 tag "ALL"
 publishDir "${OUTDIR}/Annotation/VEP", mode: 'copy'
 
input:
   file(vcf_file) from inputVep

 output:
   file(annotated_vcf) into outputVep

 when:
 	params.effect_prediction == true

 script:
   annotated_vcf = run_name + ".annotation.vep.vcf"

   """
      vep --offline --cache --dir $VEP_CACHE --fork ${task.cpus} \
 	--assembly GRCh37 -i $vcf_file -o $annotated_vcf --allele_number --canonical \
	--force_overwrite --vcf --no_progress \
	--merged \
	--pubmed \
	--plugin LoFtool --plugin LoF \
	--fasta ${params.vep_fasta}
   """

}

process runSplitVEPBySample {

        tag "ALL|${params.assembly}"
        publishDir "${OUTDIR}/Annotation/VEP/BySample", mode: 'copy'

        input:
        file(vcf) from outputVep

        output:
        file("*.vcf") into VepVcfBySample

        script:

        """
                for sample in `bcftools query -l ${vcf}`; do gatk SelectVariants -R $REF -V ${vcf} -sn \$sample -O \$sample'.vcf' ; done;
        """

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  if (params.email) {

            def subject = 'Diagnostic exome analysis finished.'
            def recipient = params.email

            ['mail', '-s', subject, recipient].execute() << """

            Pipeline execution summary
            ---------------------------
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
            """
  }

}


//#############################################################################################################
//#############################################################################################################
//
// FUNCTIONS
//
//#############################################################################################################
//#############################################################################################################


// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

