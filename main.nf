#!/usr/bin/env nextflow

/**
===============================
Exome Pipeline
===============================

This Pipeline performs one of two workflows to generate variant calls and effect predictions
using either the GATK processing chain or Freebayes.

### Homepage / git
git@github.com:marchoeppner/exome-seq.git
### Implementation
Implemented in Q1 2019

This pipeline is based on the updated GATK best-practices (where applicable).
 - trimming (FastP)
 - Alignment (BWA)
 - Duplicate marking (GATK)
 - recalibration 
 - variant calling
 - variant recalibration and filtering
 - variant effect prediction

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB Diagnostic Exome pipeline | version ${params.version}
===============================================================================
Usage: nextflow -c /path/to/git/nextflow.config run /path/to/git/main.nf --assembly hg19_clinical --kit Nextera --samples Samples.csv
This example will perform an exome analysis against the hg19 (with decoys) assembly, assuming that exome reads were generated with
the Nextera kit and using the GATK4 best-practice workflow. 
Required parameters:
--samples                      A sample list in CSV format (see website for formatting hints)
--assembly                     Name of the reference assembly to use
--kit			       Name of the exome kit (available options: xGen, xGen_custom, xGen_v2, Nextera, Pan_cancer)
--email 		       Email address to send reports to (enclosed in '')
Optional parameters:
--skip_multiqc		       Don't attached MultiQC report to the email. 
--vqsr 			       Whether to also run variant score recalibration (only works >= 30 samples) (default: false)
--panel 		       Gene panel to check coverage of (valid options: cardio_dilatative, cardio_hypertrophic, cardio_non_compaction, eoIBD_25kb, imm_eoIBD_full, breast_cancer)
--all_panels 		       Run all gene panels defined for this assembly (none if no panel is defined!)
--panel_intervals	       Run a custom gene panel in interval list format (must have a matching sequence dictionary!)
--run_name 		       A descriptive name for this pipeline run
--cram			       Whether to output the alignments in CRAM format (default: bam)
--deepvariant		       Enable variant calling with Google DeepVariant
--interval_padding	       For GATK, include this number of nt upstream and downstream around the exome targets (default: 10)
Expert options (usually not necessary to change!):
--fasta                        A reference genome in FASTA format (set automatically if using --assembly)
--dict                         A sequence dictionary matching --fasta (set automatically if using --assembly)
--dbsnp                        dbSNP data in VCF format (set automatically if using --assembly)
--g1k                          A SNP reference (usually 1000genomes, set automatically if using --assembly)
--mills_indels                 An INDEL reference (usually MILLS/1000genomes, set automatically if using --assembly)
--omni                         A SNP reference (usually OMNI, set automatically if using --assembly)
--hapmap                       A SNP reference (usually HAPMAP, set automatically if using --assembly)
--targets                      A interval_list target file (set automatically if using the --kit option)
--baits                        A interval_list bait file (set automatically if using the --kit option)
--bed                          A list of calling intervals to be used by Deepvariant (default: exome kit targets will be converted to bed)
--max_length                   Cut reads down to this length (optional, default 0 = no trimming)
--kill                         A list of known bad exons in the genome build that are ignored for panel coverage statistics (see documentation for details)
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

// #############
// INPUT OPTIONS
// #############

// Sample input file
inputFile = file(params.samples)

// Giving this pipeline run a name
params.run_name = false
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (params.run_name == false) {
	log.info "No run name was specified, using ${run_name} instead"
}

// This will eventually enable switching between multiple assembly versions
// Currently, only hg19 has all the required reference files available
params.assembly = "hg19"

REF = params.fasta ?: file(params.genomes[ params.assembly ].fasta)
DICT = params.dict ?: file(params.genomes[ params.assembly ].dict)
DBSNP = params.dbsnp ?: file(params.genomes[ params.assembly ].dbsnp )
G1K = params.g1k ?: file(params.genomes[ params.assembly ].g1k )
MILLS = params.mills_indels ?: file(params.genomes[ params.assembly ].mills )
OMNI = params.omni ?: file(params.genomes[ params.assembly ].omni )
HAPMAP = params.hapmap ?: file(params.genomes[ params.assembly ].hapmap )
// This is usually missing from target/bait definitions, so we add it 
MITOCHONDRION = params.mitochondrion ?: params.genomes[ params.assembly ].mitochondrion

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits
targets_to_bed = Channel.fromPath(TARGETS)

if (params.kill) {
	KILL = params.kill
} else if (params.genomes[params.assembly].kits[params.kit].kill) {
	KILL = params.genomes[params.assembly].kits[params.kit].kill
} else {
	KILL = false
}

SNP_RULES = params.snp_filter_rules
INDEL_RULES = params.indel_filter_rules

/*
PANEL COVERAGE - pick the correct panel for reporting
*/

if (params.panel && params.all_panels || params.panel && params.panel_intervals || params.all_panels && params.panel_intervals) {
	log.info "The options for panel stats are mutually exclusive! Will use the highest ranked choice (panel > panel_intervals > all panels)"
}
if (params.panel) {
	panel = params.genomes[params.assembly].panels[params.panel].intervals
	panels = Channel.fromPath(panel)
} else if (params.panel_intervals) {
	Channel.fromPath(params.panel_intervals)
	.ifEmpty { exit 1; "Could not find the specified gene panel (--panel_intervals)" }
	.set { panels }
} else if (params.all_panels) {
	panel_list = []
	panel_names = params.genomes[params.assembly].panels.keySet()
	panel_names.each {
		interval = params.genomes[params.assembly].panels[it].intervals
		panel_list << file(interval)
	}
	panels = Channel.fromList(panel_list)
} else {
	panels = Channel.empty()
}

// A single exon only covered in male samples - simple sex check
SRY_BED  = params.sry_bed ?: params.genomes[params.assembly].sry_bed
SRY_REGION = file(SRY_BED)

if (!SRY_REGION.exists() ) {
	exit 1, "Could not find the bed file for SRY check!"
}

// Annotations to use for variant recalibration
snp_recalibration_values = params.snp_recalibration_values
indel_recalbration_values = params.indel_recalbration_values

params.hard_filter = false

// Whether to produce BAM output instead of CRAM
params.cram = false
align_suffix = (params.cram == false) ? "bam" : "cram"

// Location of applications used
OUTDIR = file(params.outdir)

// DeepVariant variables
//params.fasta = params.genome ? params.genomes[ params.genome ].fasta : false
model = "wes"
params.fai = false
params.fastagz = false
params.gzfai = false
params.gzi = false 

// Available exome kits

if (TARGETS == false || BAITS == false ) {
   exit 1, "Information on enrichment kit incomplete or missing (please see the documentation for details!)"
}

// Whether to send a notification upon workflow completion
params.email = false

if(params.email == false) {
	exit 1, "You must provide an Email address to which pipeline updates are send!"
}

if (params.no_dedup) {
	println "Selected to skip duplicate marking. This is NOT recommended!"
}

// Check if the max_length argument is a number
if (params.max_length) {
	summary['maxReadLength'] = params.max_length
 	if (params.max_length instanceof Integer) {

	} else {
		exit 1, "Defined a max_length option, but it is not an integer...!"
	}
}

	
summary['runName'] = run_name
summary['Samples'] = inputFile
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Assembly'] = REF
summary['Kit'] = TARGETS
if (params.panel) {
	summary['GenePanel'] = params.panel
} else if (params.panel_intervals) {
	summary['GenePanel'] = params.panel_intervals
} else if (params.all_panels) {
	summary['GenePanel'] = "All panels"
}

if (KILL) {
        summary['KillList'] = KILL
}
if (workflow.containerEngine) {
	summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['References'] = [:]
summary['References']['DBSNP'] = DBSNP
summary['References']['G1K'] = G1K
summary['References']['MILLS'] = MILLS
summary['References']['OMNI'] = OMNI
summary['References']['HAPMAP'] = HAPMAP
summary['Filtering'] = [:]
summary['Filtering']['SNP_RULES'] = SNP_RULES
summary['Filtering']['INDEL_RULES'] = INDEL_RULES
summary['IntervallPadding'] = params.interval_padding
summary['SessionID'] = workflow.sessionId

// Header log info
log.info "========================================="
log.info "Exome-seq pipeline v${params.version}"
log.info "Nextflow Version:		$workflow.nextflow.version"
log.info "Assembly version: 		${params.assembly}"
log.info "Command Line:			$workflow.commandLine"
log.info "Run name: 			${run_name}"
if (workflow.containerEngine) {
	log.info "Container engine:		${workflow.containerEngine}"
}
log.info "========================================="

// Read sample file 
Channel.from(inputFile)
       .splitCsv(sep: ';', header: true)
       .set {  readPairsFastp }

process runFastp {

	scratch true

	input:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, center, date, fastqR1, fastqR2 from readPairsFastp

	output:
	set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, date, center, file(left),file(right) into inputBwa
	set file(html),file(json) into fastp_results

	script:
	def options = ""
	if (params.max_length != false) {
		options += " -b ${params.max_length} -B ${params.max_length}"
	}

	left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
	right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
	json = file(fastqR1).getBaseName() + ".fastp.json"
	html = file(fastqR1).getBaseName() + ".fastp.html"

	"""
		fastp $options --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
	"""
}

process runBWA {

    // publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Libraries/${libraryID}/${rgID}/BWA/", mode: 'copy'

    scratch true
	
    input:
    set indivID, sampleID, libraryID, rgID, platform_unit, platform, platform_model, run_date, center,file(left),file(right) from inputBwa
    
    output:
    set indivID, sampleID, file(outfile) into runBWAOutput
    
    script:
    outfile = sampleID + "_" + libraryID + "_" + rgID + ".aligned.bam"	
    
    """
	bwa mem -H $DICT -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${REF}\\tCN:${center}" -t ${task.cpus} ${REF} $left $right | samtools sort -O bam -m 2G -@ 4 - > $outfile
    """	
}

runBWAOutput_grouped_by_sample = runBWAOutput.groupTuple(by: [0,1])

process mergeBamFiles_bySample {

	scratch true

	input:
        set indivID, sampleID, file(aligned_bam_list) from runBWAOutput_grouped_by_sample

	output:
	set indivID,sampleID,file(merged_bam),file(merged_bam_index) into mergedBamFile_by_Sample, MergedBamSkipDedup

	script:
	merged_bam = sampleID + ".merged.bam"
	merged_bam_index = merged_bam + ".bai"

	if (aligned_bam_list.size() > 1 && aligned_bam_list.size() < 1000 ) {
		"""

		    	gatk MergeSamFiles \
                	    -I ${aligned_bam_list.join(' -I ')} \
	                    -O /dev/stdout \
			    --USE_THREADING true \
                	    --SORT_ORDER coordinate |

			gatk SetNmMdAndUqTags \
				-I /dev/stdin \
				-O $merged_bam \
				-R $REF \
				--IS_BISULFITE_SEQUENCE false

			samtools index $merged_bam

		"""
	} else {
		"""
			gatk SetNmMdAndUqTags \
		                -I ${aligned_bam_list.join(' -I ') } \
                	        -O $merged_bam \
                        	-R $REF \
	                        --IS_BISULFITE_SEQUENCE false

        	        samtools index $merged_bam
		"""

	}
}

process runMarkDuplicates {

        // publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/MarkDuplicates", mode: 'copy'

        scratch true

        input:
        set indivID, sampleID, file(merged_bam),file(merged_bam_index) from mergedBamFile_by_Sample

        output:
        set indivID, sampleID, file(outfile_bam),file(outfile_bai) into MarkDuplicatesOutput, BamForMultipleMetrics, runPrintReadsOutput_for_OxoG_Metrics, runPrintReadsOutput_for_HC_Metrics, BamForDepthOfCoverage
	set file(outfile_bam), file(outfile_bai) into BamForSexCheck
	file(outfile_md5) into MarkDuplicatesMD5
	file(outfile_metrics) into DuplicatesOutput_QC

        script:
        outfile_bam = sampleID + ".dedup.bam"
        outfile_bai = sampleID + ".dedup.bai"
	outfile_md5 = sampleID + ".dedup.bam.md5"

        outfile_metrics = sampleID + "_duplicate_metrics.txt"

	"""
        	gatk --java-options "-Xmx${task.memory.toGiga()-1}G" MarkDuplicates \
                	-I ${merged_bam} \
	                -O ${outfile_bam} \
        	        -M ${outfile_metrics} \
                        --CREATE_INDEX true \
			--ASSUME_SORT_ORDER=coordinate \
			--MAX_RECORDS_IN_RAM 50000 \
			--CREATE_MD5_FILE true \
                        --TMP_DIR tmp \
			-R ${REF}
	"""

}

// a simple sex check looking at coverage of the SRY gene
process runSexCheck {


	input:
	file(bams) from BamForSexCheck.collect()

	output:
	file(sex_check_yaml) into SexChecKYaml

	script:
	sex_check_yaml = "sex_check_mqc.yaml"

	"""
		parse_sry_coverage.pl --fasta $REF --region $SRY_REGION > $sex_check_yaml
	"""	
}

// If we don't want to use the deduped BAM file, allow skipping it. 
// We still run dedup for the quality stats
if (params.no_dedup) {
	MergedBamSkipDedup
        .into { BamToBSQR; BamToDV }
} else {
	MarkDuplicatesOutput
        .into { BamToBSQR; BamToDV }
}

// ************************
// Calling with DeepVariant
// ************************
if (params.deepvariant) {
	//setup fasta channels
	(fastaToIndexCh, fastaToGzCh, fastaToGzFaiCh, fastaToGziCh) = Channel.fromPath(REF).into(4)

	if(params.bed){
	bedToExamples = Channel
	    .fromPath(params.bed)
	    .ifEmpty { exit 1, "please specify --bed option (--bed bedfile)"}
	} else {
		process list_to_bed {

                	input:
	                file(targets) from targets_to_bed

        	        output:
                	file(bed) into bedToExamples

	                script:
        	        bed = targets.getBaseName() + ".bed"

                	"""
                        picard IntervalListToBed I=$targets O=$bed
			"""
		}
	}

	if(params.fai){
	faiToExamples = Channel
	    .fromPath(params.fai)
	    .ifEmpty{exit 1, "Fai file not found: ${params.fai}"}
	}

	if(params.fastagz){
	fastaGz = Channel
	    .fromPath(params.fastagz)
	    .ifEmpty{exit 1, "Fastagz file not found: ${params.fastagz}"}
	    .into {fastaGzToExamples; fastaGzToVariants }
	}

	if(params.gzfai){
	gzFai = Channel
	    .fromPath(params.gzfai)
	    .ifEmpty{exit 1, "gzfai file not found: ${params.gzfai}"}
	    .into{gzFaiToExamples; gzFaiToVariants }
	}

	if(params.gzi){
	gzi = Channel
	    .fromPath(params.gzi)
	    .ifEmpty{exit 1, "gzi file not found: ${params.gzi}"}
	    .into {gziToExamples; gziToVariants}
	}

	if(!params.fai) {
	  process preprocess_fai {
	      tag "${fasta}.fai"
	      //publishDir "${params.outdir}/DV/"

	      input:
	      file(fasta) from fastaToIndexCh

	      output:
	      file("${fasta}.fai") into faiToExamples

	      script:
	      """
	      samtools faidx $fasta
	      """
	  }
	}

	if(!params.fastagz) {
	  process preprocess_fastagz {
	      tag "${fasta}.gz"
	      //publishDir "${params.outdir}/DV/"

	      input:
	      file(fasta) from fastaToGzCh

	      output:
	      file("*.gz") into (tmpFastaGzCh, fastaGzToExamples, fastaGzToVariants)

	      script:
	      """
	      bgzip -c ${fasta} > ${fasta}.gz
	      """
	  }
	}

	if(!params.gzfai) {
	  process preprocess_gzfai {
	    tag "${fasta}.gz.fai"
	    //publishDir "${params.outdir}/DV/"

	    input:
	    file(fasta) from fastaToGzFaiCh
	    file(fastagz) from tmpFastaGzCh

	    output:
	    file("*.gz.fai") into (gzFaiToExamples, gzFaiToVariants)

	    script:
	    """
	    samtools faidx $fastagz
	    """
	  }
	}

	if(!params.gzi){
	  process preprocess_gzi {
	    tag "${fasta}.gz.gzi"
	    //publishDir "${params.outdir}/DV/"

	    input:
	    file(fasta) from fastaToGziCh

	    output:
	    file("*.gz.gzi") into (gziToExamples, gziToVariants)

	    script:
	    """
	    bgzip -c -i ${fasta} > ${fasta}.gz
	    """
	  }
	}

	/********************************************************************
	  process make_examples
	  Getting bam files and converting them to images ( named examples )
	********************************************************************/

	process DV_make_examples{

	  label 'deepvariant'

	  publishDir "${params.outdir}/make_examples", mode: 'copy',
	  saveAs: {filename -> "logs/log"}

	  input:
	  file fai from faiToExamples.collect()
	  file fastagz from fastaGzToExamples.collect()
	  file gzfai from gzFaiToExamples.collect()
	  file gzi from gziToExamples.collect()
	  file bed from bedToExamples.collect()
	  set val(indivID),val(sampleID),file(bam), file(bai) from BamToDV

	  output:
	  set val(indivID),val(sampleID),file("${bam}"),file('*_shardedExamples') into examples

	  script:
	  """
	  unset TMPDIR
	  mkdir logs
	  mkdir ${bam.baseName}_shardedExamples
	  dv_make_examples.py \
	  --cores ${task.cpus} \
	  --sample ${bam} \
	  --ref ${fastagz} \
	  --reads ${bam} \
	  --regions ${bed} \
	  --logdir logs \
	  --examples ${bam.baseName}_shardedExamples
	  """
	}

	/********************************************************************
	  process call_variants
	  Doing the variant calling based on the ML trained model.
	********************************************************************/

	process DV_call_variants{

          label 'deepvariant'

	  input:
	  set val(indivID),val(sampleID),file(bam),file(shardedExamples) from examples

	  output:
	  set val(indivID),val(sampleID),file(bam),file('*_call_variants_output.tfrecord') into called_variants

	  script:
	  """
	  dv_call_variants.py \
	    --cores ${task.cpus} \
	    --sample ${bam} \
	    --outfile ${bam.baseName}_call_variants_output.tfrecord \
	    --examples $shardedExamples \
	    --model ${model}
	  """
	}

	/********************************************************************
	  process postprocess_variants
	  Trasforming the variant calling output (tfrecord file) into a standard vcf file.
	********************************************************************/

	process DV_postprocess_variants{

          label 'deepvariant'

	  publishDir "${params.outdir}/${indivID}/${sampleID}/DeepVariant", mode: 'copy'

	  input:
	  file fastagz from fastaGzToVariants.collect()
	  file gzfai from gzFaiToVariants.collect()
	  file gzi from gziToVariants.collect()
	  set val(indivID),val(sampleID),file(bam),file('call_variants_output.tfrecord') from called_variants

	  output:
	   set val("${bam}"),file("${bam}.vcf") into postout

	  script:
	  """
	  dv_postprocess_variants.py \
	  --ref ${fastagz} \
	  --infile call_variants_output.tfrecord \
	  --outfile "${bam}.vcf"
	  """
	}

} // end Deepvariant

// ------------------------------------------------------------------------------------------------------------	//
// Perform base quality score recalibration (BQSR) including
// 1) Generate a recalibration table
// 2) Generate a new table after applying recalibration
// 3) Compare differences between recalibration tables
// 4) Apply recalibration
//
// ------------------------------------------------------------------------------------------------------------

process runBaseRecalibrator {

	// publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/BaseRecalibrator/", mode: 'copy'
    
	input:
	set indivID, sampleID, dedup_bam, dedup_bai from BamToBSQR
   
	output:
	set indivID, sampleID, dedup_bam, file(recal_table) into runBaseRecalibratorOutput

	script:
	recal_table = sampleID + "_recal_table.txt" 

	"""
		gatk --java-options "-Xmx${task.memory.toGiga()}G" BaseRecalibrator \
		--reference ${REF} \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--use-original-qualities \
		--input ${dedup_bam} \
		--known-sites ${MILLS} \
		--known-sites ${DBSNP} \
		--output ${recal_table}
	"""
}

process runApplyBQSR {

	publishDir "${OUTDIR}/${indivID}/${sampleID}/", mode: 'copy'

	scratch true
	    
	input:
	set indivID, sampleID, realign_bam, recal_table from runBaseRecalibratorOutput 

	output:
	set indivID, sampleID, file(outfile_bam), file("*.bai") into runPrintReadsOutput_for_Multiple_Metrics,inputHCSample,inputCollectReadCounts,inputPanelCoverage
	set indivID, sampleID, realign_bam, recal_table into runPrintReadsOutput_for_PostRecal
	set indivID, outfile_md5 into BamMD5
            
	script:

	def leading = ""
	if (params.run_name) {
		leading = "${run_name}."
	}
	outfile_bam = leading + sampleID + ".clean.${align_suffix}"
	outfile_bai = leading + sampleID + ".clean.${align_suffix}.bai"
	outfile_md5 = leading + sampleID + ".clean.${align_suffix}.md5"
           
    	"""
        	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyBQSR \
                --reference ${REF} \
		--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                --input ${realign_bam} \
		--use-original-qualities \
		-OBI true \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
                -bqsr ${recal_table} \
                --output ${outfile_bam} \
                -OBM true \
    	"""
}    

// Call variants on a per-sample basis

process runHCSample {

	scratch true 

	publishDir "${OUTDIR}/${indivID}/${sampleID}/Variants/HaplotypeCaller" , mode: 'copy'

	input: 
	set indivID,sampleID,file(bam),file(bai) from inputHCSample

	output:
	file(vcf) into outputHCSample
        file(vcf_index) into outputHCSampleIndex

	script:
 
	vcf = sampleID + ".raw_variants.g.vcf.gz"
	vcf_index = vcf + ".tbi"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" HaplotypeCaller \
		-R $REF \
		-I ${bam} \
		-L $TARGETS \
		-L $MITOCHONDRION \
		--dbsnp $DBSNP \
		-ip ${params.interval_padding} \
		--emit-ref-confidence GVCF \
		-OVI true \
    		--output $vcf \
		--native-pair-hmm-threads 4 \
		-G AS_StandardAnnotation  -G StandardAnnotation -G StandardHCAnnotation 
  	"""
}

// Import individual vcf files into a GenomicsDB database on a per chromosome basis
// From here on all samples are in the same file
process runGenomicsDBImport  {

        publishDir "${OUTDIR}/Variants/JointGenotypes/", mode: 'copy'

	scratch true

	input:
        file(vcf_list) from outputHCSample.collect()
	file(index_list) from outputHCSampleIndex.collect()

	output:
        set file(merged_vcf),file(merged_vcf_index) into inputGenotypeGVCFs

	script:
	merged_vcf = "merged.g.vcf.gz"
	merged_vcf_index = merged_vcf + ".tbi"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" CombineGVCFs  \
		--variant ${vcf_list.join(" --variant ")} \
		--reference $REF \
		--intervals $TARGETS \
		--intervals $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--OVI true \
		--output $merged_vcf \
		-D $DBSNP
	"""

}

// Perform genotyping on a per chromosome basis

process runGenotypeGVCFs {

	scratch true
  
	publishDir "${OUTDIR}/Variants/JointGenotypes", mode: 'copy'
  
	input:
	set file(merged_vcf), file(merged_vcf_index) from inputGenotypeGVCFs
  
	output:
	set file(gvcf), file(gvcf_index) into (inputHardFilterSNP, inputRecalSNP, inputHardFilterIndel, inputRecalIndel )
  
	script:
  
	gvcf = "genotypes.vcf.gz"
	gvcf_index = gvcf + ".tbi"
  
	"""
 	gatk --java-options "-Xmx${task.memory.toGiga()}G" GenotypeGVCFs \
		--reference $REF \
		--dbsnp $DBSNP \
		-L $TARGETS \
		-L $MITOCHONDRION \
		-ip ${params.interval_padding} \
		--only-output-calls-starting-in-intervals \
		-V $merged_vcf \
              	--output $gvcf \
                -G StandardAnnotation -G AS_StandardAnnotation \
		-OVI true
	"""
}

////////////////////////
// Hard filtering
////////////////////////

process runHardFilterSNP {
		
	publishDir "${OUTDIR}/Variants/HardFilter/Preprocess", mode: 'copy'

	input:
	set file(vcf),file(vcf_index) from inputHardFilterSNP

	output:
	set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterSNP

	script:
	vcf_filtered = "genotypes.merged.snps.filtered.vcf.gz"
	vcf_filtered_index = vcf_filtered + ".tbi"

	"""
		gatk SelectVariants \
			-R $REF \
			-V $vcf \
			--select-type-to-include SNP \
			--output genotypes.merged.snps.vcf.gz \
			-OVI true
				
		gatk VariantFiltration \
			-R $REF \
			-V genotypes.merged.snps.vcf.gz \
			--output $vcf_filtered \
			--filter-expression "${SNP_RULES}" \
			--filter-name "hard_snp_filter" \
			-OVI true
	"""

}

process runHardFilterIndel {

        publishDir "${OUTDIR}/Variants/HardFilter/Preprocess", mode: 'copy'
        
        input:
        set file(vcf),file(vcf_index) from inputHardFilterIndel

        output:
        set file(vcf_filtered),file(vcf_filtered_index) into outputHardFilterIndel

        script:
        vcf_filtered = "genotypes.merged.indels.filtered.vcf.gz"
        vcf_filtered_index = vcf_filtered + ".tbi"

        """
 	       gatk SelectVariants \
               -R $REF \
               -V $vcf \
               --select-type-to-include INDEL \
               --output genotypes.merged.indels.vcf.gz \
               -OVI true

              gatk VariantFiltration \
              -R $REF \
              -V genotypes.merged.indels.vcf.gz \
              --output $vcf_filtered \
	      --filter-expression "${INDEL_RULES}" \
              --filter-name "hard_indel_filter" \
	      -OVI true
        """
}

process runCombineHardVariants {

        publishDir "${OUTDIR}/Variants/HardFilter/Final", mode: 'copy'

        input:
        set file(indel),file(indel_index) from outputHardFilterIndel
	set file(snp),file(snp_index) from outputHardFilterSNP

        output:
        set file(merged_file),file(merged_file_index) into inputfilterPassVariants

        script:
        merged_file = "${run_name}.merged_callset.hard_filter.vcf.gz"
        merged_file_index = merged_file + ".tbi"

        """
		gatk --java-options "-Xmx${task.memory.toGiga()}G"  MergeVcfs \
		-I $indel \
		-I $snp \
		-O $merged_file \
		
		gatk IndexFeatureFile -I $merged_file

        """
}

process runSplitHardVariantsBySample {

        publishDir "${OUTDIR}/Variants/HardFilter/Final/BySample", mode: 'copy'

        input:
        set file(vcf_clean),file(vcf_clean_index) from inputfilterPassVariants

        output:
        file("*.vcf.gz*") into HardVcfBySample

        script:

        """
                for sample in `bcftools query -l $vcf_clean`; do gatk SelectVariants -R $REF -V $vcf_clean --exclude-non-variants --remove-unused-alternates -sn \$sample -O \$sample'.vcf.gz' ; done;
        """

}

/////////////////////////
// Variant recalibration 
/////////////////////////

process runRecalibrationModeSNP {

	publishDir "${OUTDIR}/Variants/VSQR/Recal"
	
	input:
	set file(vcf),file(vcf_index) from inputRecalSNP

	output:
	set file(recal_file),file(tranches) into inputRecalSNPApply

	when:
	params.vqsr == true

	script:
	recal_file = "genotypes.merged.snps.recal"
	tranches = "genotypes.merged.snps.tranches"

	"""

		gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
		-R $REF \
		-V $vcf \
	       	-O $recal_file \
       	        --tranches-file $tranches \
		-an ${snp_recalibration_values.join(' -an ')} \
	        -mode SNP \
		-OVI true \
		--resource hapmap,known=false,training=true,truth=true,prior=15.0:$HAPMAP \
		--resource omni,known=false,training=true,truth=true,prior=12.0:$OMNI \
		--resource 1000G,known=false,training=true,truth=false,prior=10.0:$G1K \
		--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
		-tranche ${params.snp_recalibration_tranche_values.join(' -tranche ')} \
		--max-gaussians 4
	"""
}

process runRecalibrationModeIndel {
	
	publishDir "${OUTDIR}/Variants/VSQR/Recal"

  	input:
	set file(vcf),file(vcf_index) from inputRecalIndel

  	output:
	set file(recal_file),file(tranches),file(vcf),file(vcf_index) into inputRecalIndelApply

	when:
	params.vqsr == true

	script:
  	recal_file = "genotypes.merged.indel.recal"
	tranches = "genotypes.merged.indel.tranches"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantRecalibrator \
        	        -R $REF \
	        	-V $vcf \
	               	-O $recal_file \
        	        --tranches-file $tranches \
	                -an ${indel_recalibration_values.join(' -an ')} \
	                -mode INDEL \
			-OVI true \
	        	--resource mills,known=false,training=true,truth=true,prior=15.0:$MILLS \
	               	--resource dbsnp,known=true,training=false,truth=false,prior=2.0:$DBSNP \
			-tranche ${params.indel_recalibration_tranche_values.join(' -tranche ')} \
			--max-gaussians 3
  	"""
}

process runRecalIndelApply {

        publishDir "${OUTDIR}/Variants/VSQR/Recal"

        input:
        set file(recal_file),file(tranches),file(gvcf),file(gvcf_index) from inputRecalIndelApply

        output:
        set file(vcf_indel),file(vcf_indel_index) into outputRecalIndelApply

        script:

        vcf_indel = "genotypes.recal_Indel.vcf.gz"
        vcf_indel_index = vcf_indel + ".tbi"

        """
        	gatk IndexFeatureFile -F $recal_file
                        gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
                        -R $REF \
                        -V $gvcf \
                        --recal-file $recal_file \
                        --tranches-file $tranches \
                        -mode INDEL \
                        --ts-filter-level ${params.indel_filter_level} \
                        -OVI true \
                         -O $vcf_indel
        """
}

process runRecalSNPApply {
	
	publishDir "${OUTDIR}/Variants/VSQR/Filtered"
	
	input:
	set file(vcf),file(index) from outputRecalIndelApply
	set file(recal_file),file(tranches) from inputRecalSNPApply

	output:
	set file(vcf_snp),file(vcf_snp_index) into outputRecalSNPApply

	script:
 
	vcf_snp = "genotypes.recal_Indel.recal_SNP.vcf.gz"
	vcf_snp_index = vcf_snp + ".tbi"

	"""
	gatk IndexFeatureFile -F $recal_file
	gatk --java-options "-Xmx${task.memory.toGiga()}G" ApplyVQSR \
		-R $REF \
		-V $vcf \
	        --recal-file $recal_file \
       	       	--tranches-file $tranches \
		-mode SNP \
		--ts-filter-level ${params.snp_filter_level} \
		-O $vcf_snp \
		-OVI true	
	"""
}

process runVariantFiltrationIndel {

	publishDir "${OUTDIR}/Variants/VSQR/Filtered"

  	input:
	set file(vcf),file(vcf_index) from outputRecalIndelApply

  	output:
  	set file(filtered_gvcf),file(filtered_gvcf_index) into inputSelectVariants

  	script:

  	filtered_gvcf = "genotypes.filtered.final.vcf.gz"
	filtered_gvcf_index = filtered_gvcf + ".tbi"

	"""
	gatk --java-options "-Xmx${task.memory.toGiga()}G" VariantFiltration \
       	       -R $REF \
               	-V $vcf \
		--filter-expression "QD < 2.0" \
		--filter-name "QDFilter" \
               	--output $filtered_gvcf \
		-OVI true
  	"""
}

process runSelectVariants {

	publishDir "${OUTDIR}/Variants/VSQR/Final", mode: 'copy'

	input:
	set file(vcf),file(vcf_index) from inputSelectVariants

	output:
	set file(vcf_clean),file(vcf_clean_index) into inputVep, inputSplitSample

	script:
	vcf_clean = "${run_name}.variants.merged.filtered.controls_removed.vcf.gz"
	vcf_clean_index = vcf_clean + ".tbi"

	"""
		gatk SelectVariants \
		-V $vcf \
		-R $REF \
		-O $vcf_clean \
		--remove-unused-alternates true \
		-OVI true \
		--exclude-non-variants true \
	"""

}

process runSplitBySample {

        publishDir "${OUTDIR}/Variants/VSQR/Final/BySample", mode: 'copy'

	input:
	set file(vcf_clean),file(vcf_clean_index) from inputSplitSample

	output: 
	file("*.vcf.gz*") into VcfBySample

	script: 

	"""
		for sample in `bcftools query -l $vcf_clean`; do gatk SelectVariants -R $REF -V $vcf_clean --exclude-non-variants --remove-unused-alternates -sn \$sample -O \$sample'.vcf.gz' -OVI true ; done;
	"""

}

// *********************
// Compute statistics for fastQ files, libraries and samples
// *********************

process runCollectMultipleMetrics {
	publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'
 
	scratch true
	    
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

    publishDir "${OUTDIR}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

    input:
    set indivID, sampleID, file(bam), file(bai) from runPrintReadsOutput_for_HC_Metrics

    output:
    file(outfile) into HybridCaptureMetricsOutput mode flatten
    file(outfile_per_target) into HsMetricsPerTarget

    script:
    outfile = sampleID + ".hybrid_selection_metrics.txt"
    outfile_per_target = sampleID + ".hybrid_selection_per_target_metrics.txt"

    """
        picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
		PER_TARGET_COVERAGE=${outfile_per_target} \
                TARGET_INTERVALS=${TARGETS} \
                BAIT_INTERVALS=${BAITS} \
                REFERENCE_SEQUENCE=${REF} \
                TMP_DIR=tmp
        """
}

process runOxoGMetrics {

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

// this is not finished yet, need to create a proper yaml file
process get_software_versions {

    publishDir "${OUTDIR}/Summary/versions", mode: 'copy'
    output:
    file("v*.txt") into software_versions
    file(yaml_file) into (software_versions_yaml_fastqc, software_versions_yaml_lib, software_versions_yaml_sample)

    script:
    yaml_file = "software_versions_mqc.yaml"

    """
    echo $workflow.manifest.version &> v_ikmb_exoseq.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    fastp -v &> v_fastp.txt
    gatk --version &> v_gatk.txt
    picard MarkDuplicates -h &> /dev/stdout | grep "Version" > v_picard.txt  || true
    samtools --version &> v_samtools.txt
    multiqc --version &> v_multiqc.txt
    parse_versions.pl >  $yaml_file
    """
}

process runMultiqcFastq {

    publishDir "${OUTDIR}/Summary/Fastq", mode: 'copy'

    when:
    !params.skip_multiqc	    

    input:
    file('*') from fastp_results.flatten().toList()
    file('*') from software_versions_yaml_fastqc.collect()
    
    output:
    file("fastp_multiqc*") into runMultiQCFastqOutput
    	
    script:

    """
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    cp $params.logo . 
    multiqc -n fastp_multiqc *.json *.html *.yaml
    """
}

process runMultiqcLibrary {

    publishDir "${OUTDIR}/Summary/Library", mode: 'copy'

    when:
    !params.skip_multiqc
	    
    input:
    file('*') from DuplicatesOutput_QC.flatten().toList()
    file('*') from software_versions_yaml_lib.collect()

    output:
    file("library_multiqc*") into runMultiQCLibraryOutput
    	
    script:

    """
    cp $params.logo .
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n library_multiqc *.txt *.yaml
    """
}

process runMultiqcSample {

    publishDir "${OUTDIR}/Summary/Sample", mode: 'copy'

    when:
    !params.skip_multiqc
	    
    input:
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file('*') from runOxoGMetricsOutput.flatten().toList()
    file('*') from software_versions_yaml_sample.collect() 
    file('*') from SexChecKYaml.collect()
        
    output:
    file("sample_multiqc.html") into multiqc_report
    file("*.yaml")
	
    script:

    def subject = 'Diagnostic exome analysis quality report'
    def recipient = params.email

    """
    cp $params.logo . 
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n sample_multiqc *

    """
}

panel_coverage_data = inputPanelCoverage.combine(panels)

process runPanelCoverage {

	publishDir "${OUTDIR}//Summary/Panel/PanelCoverage", mode: "copy"

        input:
        set indivID,sampleID,file(bam),file(bai),file(panel) from panel_coverage_data

        output:
        set val(panel_name),file(coverage) into outputPanelCoverage
	set indivID,sampleID,file(target_coverage_xls) into outputPanelTargetCoverage
	file(target_coverage)

        script:
        panel_name = panel.getSimpleName()
        coverage = indivID + "_" + sampleID + "." +  panel_name  + ".hs_metrics.txt"
	target_coverage = indivID + "_" + sampleID + "." +  panel_name  + ".per_target.hs_metrics.txt"
	target_coverage_xls = indivID + "_" + sampleID + "." + panel_name + ".per_target.hs_metrics_mqc.xlsx"

	// optionally support a kill list of known bad exons
	def options = ""
	if (KILL) {
		options = "--ban ${KILL}"
	}

        // do something here - get coverage and build an XLS sheet
	// First we identify which analysed exons are actually part of the exome kit target definition. 
        """

		picard -Xmx${task.memory.toGiga()}G IntervalListTools \
			INPUT=$panel \
			SECOND_INPUT=$TARGETS \
			ACTION=SUBTRACT \
			OUTPUT=overlaps.interval_list

                picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                        INPUT=${bam} \
                        OUTPUT=${coverage} \
                        TARGET_INTERVALS=${panel} \
                        BAIT_INTERVALS=${panel} \
			CLIP_OVERLAPPING_READS=false \
                        REFERENCE_SEQUENCE=${REF} \
                        TMP_DIR=tmp \
			PER_TARGET_COVERAGE=$target_coverage

		target_coverage2xls.pl $options --infile $target_coverage --min_cov $params.panel_coverage --skip overlaps.interval_list --outfile $target_coverage_xls

        """
}

grouped_panels = outputPanelCoverage.groupTuple()

process runMultiqcPanel {

	publishDir "${OUTDIR}/Summary/Panel", mode: "copy"

	input:
	set val(panel_name),file('*') from grouped_panels

	output:
	file("${panel_name}_multiqc.html") into panel_qc_report

	script:

	"""
		cp $params.logo . 
		cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
		multiqc -n ${panel_name}_multiqc *
	"""
}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="

  def email_fields = [:]
  email_fields['version'] = workflow.manifest.version
  email_fields['session'] = workflow.sessionId
  email_fields['runName'] = run_name
  email_fields['Samples'] = params.samples
  email_fields['success'] = workflow.success
  email_fields['dateStarted'] = workflow.start
  email_fields['dateComplete'] = workflow.complete
  email_fields['duration'] = workflow.duration
  email_fields['exitStatus'] = workflow.exitStatus
  email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
  email_fields['errorReport'] = (workflow.errorReport ?: 'None')
  email_fields['commandLine'] = workflow.commandLine
  email_fields['projectDir'] = workflow.projectDir
  email_fields['script_file'] = workflow.scriptFile
  email_fields['launchDir'] = workflow.launchDir
  email_fields['user'] = workflow.userName
  email_fields['Pipeline script hash ID'] = workflow.scriptId
  email_fields['kit'] = TARGETS
  email_fields['assembly'] = REF
  email_fields['manifest'] = workflow.manifest
  email_fields['summary'] = summary

  email_info = ""
  for (s in email_fields) {
	email_info += "\n${s.key}: ${s.value}"
  }

  def output_d = new File( "${params.outdir}/pipeline_info/" )
  if( !output_d.exists() ) {
      output_d.mkdirs()
  }

  def output_tf = new File( output_d, "pipeline_report.txt" )
  output_tf.withWriter { w -> w << email_info }	

 // make txt template
  def engine = new groovy.text.GStringTemplateEngine()

  def tf = new File("$baseDir/assets/email_template.txt")
  def txt_template = engine.createTemplate(tf).make(email_fields)
  def email_txt = txt_template.toString()

  // make email template
  def hf = new File("$baseDir/assets/email_template.html")
  def html_template = engine.createTemplate(hf).make(email_fields)
  def email_html = html_template.toString()
  
  def subject = "Diagnostic exome analysis finished ($run_name)."

  if (params.email) {

  	def mqc_report = null
  	try {
        	if (workflow.success && !params.skip_multiqc) {
            		mqc_report = multiqc_report.getVal()
            		if (mqc_report.getClass() == ArrayList){
                		log.warn "[IKMB ExoSeq] Found multiple reports from process 'multiqc', will use only one"
                		mqc_report = mqc_report[0]
                	}
        	}
    	} catch (all) {
        	log.warn "[IKMB ExoSeq] Could not attach MultiQC report to summary email"
  	}

	def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
	def sf = new File("$baseDir/assets/sendmail_template.txt")	
    	def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    	def sendmail_html = sendmail_template.toString()

	try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
        }

  }

}
