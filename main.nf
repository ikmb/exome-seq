#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
Exome Pipeline
===============================

This Pipeline performs variant calling
using Google DeepVariant

### Homepage / git
git@github.com:ikmb/exome-seq.git
### Implementation
Re-Implemented in Q3 2020

This pipeline is based on the DeepVariant best-practices (where applicable).
 - trimming (FastP)
 - Alignment (BWA)
 - Duplicate marking (Samtools)
 - Variant calling per-sample (Deepvariant)
 - Variant calling across samples (GLNexus)

Author: Marc P. Hoeppner, m.hoeppner@ikmb.uni-kiel.de

**/

// Pipeline version

def summary = [:]

// #############
// INPUT OPTIONS
// #############

// Set Channels
params.fasta = file(params.genomes[ params.assembly ].fasta, checkIfExists: true)
params.fasta_fai = file(params.genomes[ params.assembly ].fai, checkIfExists: true)
params.fasta_gz = file(params.genomes[ params.assembly ].fastagz, checkIfExists: true)
params.fasta_gzfai = file(params.genomes[ params.assembly ].gzfai, checkIfExists: true)
params.fasta_gzi = file(params.genomes[ params.assembly ].gzi, checkIfExists: true)
params.dict = file(params.genomes[ params.assembly ].dict, checkIfExists: true)
params.dbsnp = file(params.genomes[ params.assembly ].dbsnp, checkIfExists: true)
params.csq_gtf = params.genomes[params.assembly].gtf
params.omni = file( params.genomes[ params.assembly ].omni, checkIfExists: true)
params.hapmap = file(params.genomes[ params.assembly ].hapmap, checkIfExists: true)
params.g1k = file(params.genomes[ params.assembly].g1k, checkIfExists: true)
params.mills = file(params.genomes[ params.assembly ].mills, checkIfExists: true)
params.axiom = file(params.genomes[ params.assembly ].axiom, checkIfExists: true)

params.snps = [ params.hapmap,  params.omni,  params.dbsnp, params.g1k ]
params.indels = [ params.mills,  params.axiom ]

if (params.amplicon_bed) { ch_amplicon_bed = Channel.fromPath(file(params.amplicon_bed, checkIfExists: true)) } else { ch_amplicon_bed = Channel.from([]) }


ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)

deepvariant_ref = Channel.from( [ file(params.fasta_fai),file(params.fasta_gz),file(params.fasta_gzfai),file(params.fasta_gzi)] )

if (params.bwa2) {
        params.bwa_index = file(params.genomes[ params.assembly ].bwa2_index, checkIfExists: true)
} else {
        params.bwa_index = file(params.genomes[ params.assembly ].fasta, checkIfExists: true)
}

//
// Targets and Baits
//
TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

if (TARGETS==BAITS) {
        exit 1, "Target and bait files must not have the same name to avoid file collisions!"
}

targets = Channel.from(file(TARGETS, checkIfExists: true))
baits = Channel.from(file(BAITS, checkIfExists: true))

//
// Kill List
//
if (params.kill) {
        params.kill_list = params.kill
} else if (params.kit && params.genomes[params.assembly].kits[params.kit].kill) {
        params.kill_list = params.genomes[params.assembly].kits[params.kit].kill
} else {
        params.kill_list = false
}

/*
PANEL COVERAGE - pick the correct panel for reporting
*/

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
params.sry_region  = params.sry_bed ?: params.genomes[params.assembly].sry_bed

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

//
// Input validation
//
WorkflowMain.initialise(workflow, params, log)
WorkflowExomes.initialise( params, log)

// Expansion hunter references
if ('expansionhunter' in tools) {
        ecatalog = file(params.genomes[params.assembly].expansion_catalog, checkIfExists: true )
        Channel.fromPath(ecatalog)
        .ifEmpty { exit 1, "Could not find a matching ExpansionHunter catalog for this assembly" }
        .set { expansion_catalog }
} else {
        expansion_catalog = Channel.empty()
}

if ('cnvkit' in tools) {
        params.cnv_ref = params.cnv_gz ?: file(params.genomes[params.assembly ].kits[params.kit].cnv_ref)
} else {
        params.cnv_ref = Channel.empty()
}

//
// Summary of all options
//
summary['runName'] = params.run_name
summary['Samples'] = params.samples
summary['Current home'] = "$HOME"
summary['Current user'] = "$USER"
summary['Current path'] = "$PWD"
summary['Assembly'] = params.fasta
if (params.bwa2) {
        summary['ALIGNER'] = "bwa-mem2"
} else {
        summary['ALIGNER'] = "bwa"
}
summary['Kit'] = TARGETS
if (params.panel) {
        summary['GenePanel'] = params.panel
} else if (params.panel_intervals) {
        summary['GenePanel'] = params.panel_intervals
} else if (params.all_panels) {
        summary['GenePanel'] = "All panels"
}
summary['JointCalling'] = params.joint_calling
if ('expansionhunter' in tools) {
	summary['ExpansionAnalysis'] = true
}
summary['CommandLine'] = workflow.commandLine
if (params.kill_list) {
        summary['KillList'] = params.kill_list
}
if (workflow.containerEngine) {
        summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['References'] = [:]
summary['References']['DBSNP'] = params.dbsnp
if (params.effects) {
        summary['References']['dbNSFP'] = params.dbnsfp_db
        summary['References']['dbSCSNV'] = params.dbscsnv_db
        summary['References']['CADD_SNPs'] = params.cadd_snps
        summary['References']['CADD_Indels'] = params.cadd_indels
}
if ('cnvkit' in tools) {
	summary['CNVkit'] = [:]
	summary['Reference'] = params.cnv_ref
}
summary['IntervallPadding'] = params.interval_padding
summary['SessionID'] = workflow.sessionId
if (params.amplicon_bed) {
	summary['AmpliconBedFile'] = params.amplicon_bed
}
// Read sample file
ch_samplesheet = file(params.samples, checkIfExists: true)

//
// Modules and workflows to include
//
include { CONVERT_BED } from "./workflows/bed/main.nf"
include { TRIM_AND_ALIGN } from "./workflows/align/main.nf"
include { DV_VARIANT_CALLING } from "./workflows/deepvariant/main.nf"
include { GATK_VARIANT_CALLING } from "./workflows/gatk/main.nf"
include { STRELKA_VARIANT_CALLING ; STRELKA_MULTI_CALLING } from "./workflows/strelka/main.nf"
include { MERGE_GVCFS } from "./modules/deepvariant/main.nf"
include { MANTA } from "./modules/manta/main.nf"
include { PICARD_METRICS } from "./workflows/picard/main.nf"
include { EXPANSIONS } from "./workflows/expansions/main.nf"
include { VEP; HAPLOSAURUS } from "./modules/vep/main.nf"
include { MULTIQC as  multiqc_fastq ; MULTIQC as multiqc_library ; MULTIQC as multiqc_sample } from "./modules/multiqc/main.nf"
include { MERGE_VCF ; VCF_ADD_DBSNP; VCF_STATS ; VCF_INDEX } from "./modules/vcf/main.nf"
include { PANEL_QC } from "./workflows/panels/main.nf"
include { CNVKIT } from "./workflows/cnvkit/main.nf"
include { CSQ; CONCAT } from "./modules/bcftools/main.nf" 
include { SEX_CHECK} from "./modules/qc/main.nf"

def multiqc_report = []

workflow {

	main:

                ch_vcfs = Channel.empty()
		ch_phased_vcfs = Channel.empty()

		// create calling regions
		CONVERT_BED(targets)
		padded_bed = CONVERT_BED.out.bed
		bedgz = CONVERT_BED.out.bed_gz

		// align reads against genome
		TRIM_AND_ALIGN(ch_samplesheet,ch_amplicon_bed)
		bam = TRIM_AND_ALIGN.out.bam
		trim_report = TRIM_AND_ALIGN.out.qc
		dedup_report = TRIM_AND_ALIGN.out.dedup_report
		sample_names = TRIM_AND_ALIGN.out.sample_names

		// DEEPVARIANT WORKFLOW
		if ('deepvariant' in tools) {
			DV_VARIANT_CALLING(bam,padded_bed,deepvariant_ref.collect())
			dv_vcf = DV_VARIANT_CALLING.out.vcf
                	dv_merged_vcf = DV_VARIANT_CALLING.out.vcf_multi
			ch_vcfs = ch_vcfs.mix(dv_vcf,dv_merged_vcf)
			ch_phased_vcfs = ch_phased_vcfs.mix(DV_VARIANT_CALLING.out.vcf_phased_multi,DV_VARIANT_CALLING.out.vcf_phased_single)
		} else {
			dv_vcf = Channel.empty()
			dv_merged_vcf = Channel.empty()
		}

		// GATK WORKFLOW
		if ('gatk' in tools) {
			GATK_VARIANT_CALLING(
				bam,
				targets,
				TRIM_AND_ALIGN.out.metas
			)
			gatk_vcf = GATK_VARIANT_CALLING.out.vcf
			gatk_merged_vcf = GATK_VARIANT_CALLING.out.vcf_multi
			ch_vcfs = ch_vcfs.mix(GATK_VARIANT_CALLING.out.vcf_multi)
			
		} else {
			gatk_vcf = Channel.empty()
			gatk_merged_vcf = Channel.empty()
		}

		// STRELKA WORKFLOW
		if ('strelka' in tools) {
			// Call all samples together
                        if (params.joint_calling) {
				STRELKA_MULTI_CALLING(
					bam.map {m,b,i -> [ b,i ] },
					bedgz,
					TRIM_AND_ALIGN.out.metas
				)
				ch_vcfs = ch_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf,STRELKA_MULTI_CALLING.out.vcf_multi)
				ch_phased_vcfs = ch_phased_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf_phased_multi)
			// Call each sample individually and merge later
                        } else {
	        	        STRELKA_VARIANT_CALLING(bam,bedgz,sample_names)
        	        	strelka_vcf = STRELKA_VARIANT_CALLING.out.vcf
				strelka_merged_vcf = STRELKA_VARIANT_CALLING.out.vcf_multi
				ch_vcfs = ch_vcfs.mix(strelka_vcf,strelka_merged_vcf)
			}
		} else {
			strelka_vcf = Channel.empty()
			strelka_merged_vcf = Channel.empty()
		}
		
		// CNV Calling
		if ('cnvkit' in tools) {
			CNVKIT(padded_bed,bam,file(params.cnv_ref))
		}

		// SV calling with Manta
		if ('manta' in tools) {
			MANTA(bam,bedgz.collect())
			manta_vcf = MANTA.out[0].mix(MANTA.out[1],MANTA.out[2])
		} else {
			manta_vcf = Channel.empty()
		}

		// Expansions
                if ('expansionhunter' in tools) {
	 		EXPANSIONS(bam,expansion_catalog)
		}

		// QC Metrics
		PANEL_QC(bam,panels,targets)
		PICARD_METRICS(bam,targets,baits)
		bam_qc = PICARD_METRICS.out.qc_reports
		VCF_STATS(ch_vcfs)
		vcf_qc = VCF_STATS.out.stats

		// Coverage of SRY gene
		SEX_CHECK(
			bam.map { m,b,i -> tuple(b,i) }.collect()
		)

		// Effect prediction
		if (params.effects) {
			VEP(ch_vcfs)
			CSQ(ch_phased_vcfs)
			if ('haplosaurus' in tools) {
				HAPLOSAURUS(ch_phased_vcfs)
			}
		}

		// QC Reports
		multiqc_fastq("FastQ",trim_report.collect())
		multiqc_library("Library",dedup_report.collect())
		multiqc_sample("Sample",bam_qc.mix(vcf_qc,SEX_CHECK.out.yaml).collect())

                multiqc_report = multiqc_sample.out.report.toList()
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
  email_fields['assembly'] = FASTA
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
  
  def subject = "Diagnostic exome analysis finished ($params.run_name)."

  if (params.email) {

        def mqc_report = null
        try {
                if (workflow.success && !params.skip_multiqc) {
                        mqc_report = multiqc_report.getVal()
                        if (mqc_report.getClass() == ArrayList){
                                log.warn "[IKMB ExoSeq] Found multiple reports from process 'multiqc', will use only one"
                                mqc_report = mqc_report[-1]
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

