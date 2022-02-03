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

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
IKMB Diagnostic Exome pipeline | version ${params.version}
===============================================================================
Usage: nextflow run ikmb/exome-seq --assembly GRCh38 --kit xGen_v2 --samples Samples.csv
This example will perform an exome analysis against the ALT-free hg38 assembly, assuming that exome reads were generated with
the IDT xGen v2 kit and using DeepVariant with GLNexus.

Required parameters:
--samples                      A sample list in CSV format (see website for formatting hints)
--assembly                     Name of the reference assembly to use
--kit                          Name of the exome kit (available options: xGen, xGen_custom, xGen_v2, Nextera, Pan_cancer)
--email                        Email address to send reports to (enclosed in '')
Optional parameters:
--bwa2                         Use BWA2 instead of BWA1
--expansion_hunter             Run ExpansionHunter
--joint_calling                Perform joint calling of all samples (default: true)
--skip_multiqc                 Don't attached MultiQC report to the email.
--panel                        Gene panel to check coverage of (valid options: cardio_dilatative, cardio_hypertrophic, cardio_non_compaction, eoIBD_25kb, imm_eoIBD_full, breast_cancer)
--all_panels                   Run all gene panels defined for this assembly (none if no panel is defined!)
--panel_intervals              Run a custom gene panel in interval list format (must have a matching sequence dictionary!)
--run_name                     A descriptive name for this pipeline run
--cram                         Whether to output the alignments in CRAM format (default: bam)
--interval_padding             Include this number of nt upstream and downstream around the exome targets (default: 10)
--vep                          Perform variant annotation with VEP (requires substantial local configuration work!)
Expert options (usually not necessary to change!):
--fasta                        A reference genome in FASTA format (set automatically if using --assembly)
--dict                         A sequence dictionary matching --fasta (set automatically if using --assembly)
--dbsnp                        dbSNP data in VCF format (set automatically if using --assembly)
--targets                      A interval_list target file (set automatically if using the --kit option)
--baits                        A interval_list bait file (set automatically if using the --kit option)
--bed                          A list of calling intervals to be used by Deepvariant (default: exome kit targets will be converted to bed)
--min_mapq                     Minimum mapping quality to consider for general coverage analysis (default = 20)
--kill                         A list of known bad exons in the genome build that are ignored for panel coverage statistics (see documentation for details)
Output:
--outdir                       Local directory to which all output is written (default: results)
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
}

def summary = [:]

// #############
// INPUT OPTIONS
// #############

// Giving this pipeline run a name
run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

if (!params.run_name) {
        exit 1, "No run name was specified, exiting..."
}

params.fasta = file(params.genomes[ params.assembly ].fasta)
params.fasta_fai = file(params.genomes[ params.assembly ].fai)
params.fasta_gz = file(params.genomes[ params.assembly ].fastagz)
params.fasta_gzfai = file(params.genomes[ params.assembly ].gzfai)
params.fasta_gzi = file(params.genomes[ params.assembly ].gzi)
params.dict = file(params.genomes[ params.assembly ].dict)
params.dbsnp = file(params.genomes[ params.assembly ].dbsnp )
params.cnv_ref = file(params.genomes[params.assembly ].kits[params.kit].cnv_ref)

deepvariant_ref = Channel.from( [ file(params.fasta_fai),file(params.fasta_gz),file(params.fasta_gzfai),file(params.fasta_gzi)] )

if (params.cnv && !file(params.cnv_ref).exists()) {
	exit 1, "Missing cnv ref file for this kit"
}
cnv_cnn = Channel.from(params.cnv_ref)

if (params.bwa2) {
        params.bwa_index = file(params.genomes[ params.assembly ].bwa2_index)
} else {
        params.bwa_index = file(params.genomes[ params.assembly ].fasta)
}

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

if (TARGETS==BAITS) {
        exit 1, "Target and bait files must not have the same name to avoid file collisions!"
}

targets = Channel.from(file(TARGETS))
baits = Channel.from(file(BAITS))

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
params.sry_region  = params.sry_bed ?: params.genomes[params.assembly].sry_bed
SRY_REGION = file(params.sry_region)

if (!SRY_REGION.exists() ) {
        exit 1, "Could not find the bed file for SRY check!"
}

// Expansion hunter references
if (params.expansion_hunter) {
        ecatalog = file(params.genomes[params.assembly].expansion_catalog)
        Channel.fromPath(ecatalog)
        .ifEmpty { exit 1, "Could not find a matching ExpansionHunter catalog for this assembly" }
        .set { expansion_catalog }
} else {
        expansion_catalog = Channel.empty()
}

// Whether to produce BAM output instead of CRAM
params.cram = false
align_suffix = (params.cram == false) ? "bam" : "cram"

// Input validation
if (params.vep) {

        if (params.assembly != "GRCh38") {
                exit 1, "VEP is not currently set up to work with defunct assembly versions...please use GRCh38"
        }

        if (!params.vep_cache_dir || !params.vep_plugin_dir) {
                exit 1, "Missing VEP cache and/or plugin directory..."
        }

        if (params.dbnsfp_db) {
                dbNSFP_DB = file(params.dbnsfp_db)
                if (!dbNSFP_DB.exists()) {
                        exit 1, "Could not find the specified dbNSFP database..."
                }
        } else {
                exit 1, "No dbNSFP database defined for this execution profile..."
        }

        if (params.dbscsnv_db) {
                dbscSNV_DB = file(params.dbscsnv_db)
                if ( !dbscSNV_DB.exists() ) {
                        exit 1, "Could not find the specified dbscSNV database..."
                }
        } else {
                exit 1, "No dbscSNV database defined for this execution profile..."
        }

        if (params.cadd_snps && params.cadd_indels) {
                CADD_SNPS = file(params.cadd_snps)
                CADD_INDELS = file(params.cadd_indels)
                if (!CADD_SNPS.exists() || !CADD_INDELS.exists() )
                        exit 1, "Missing CADD SNPs and/or Indel references..."
                }
        else {
                exit 1, "CADD SNP and/or Indel reference files not defined for this execution profile..."
        }

}

summary['runName'] = run_name
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
summary['ExpansionAnalysis'] = params.expansion_hunter
summary['CommandLine'] = workflow.commandLine
if (params.kill_list) {
        summary['KillList'] = params.kill_list
}
if (workflow.containerEngine) {
        summary['Container'] = "$workflow.containerEngine - $workflow.container"
}
summary['References'] = [:]
summary['References']['DBSNP'] = params.dbsnp
if (params.vep) {
        summary['References']['dbNSFP'] = params.dbnsfp_db
        summary['References']['dbSCSNV'] = params.dbscsnv_db
        summary['References']['CADD_SNPs'] = params.cadd_snps
        summary['References']['CADD_Indels'] = params.cadd_indels
}
if (params.cnv) {
	summary['CNVkit'] = [:]
	summary['CNVkit']['BlackList'] = params.cnv_blacklist
	summary['CNVkit']['ExcludeRegions'] = params.cnv_exclusion
}
summary['IntervallPadding'] = params.interval_padding
summary['SessionID'] = workflow.sessionId

// Header log info
log.info "========================================="
log.info "Exome-seq pipeline v${params.version}"
log.info "Nextflow Version:             $workflow.nextflow.version"
log.info "Assembly version:             ${params.assembly}"
if (params.bwa2) {
        log.info "Read aligner:                 BWA2"
} else {
        log.info "Read aligner:                 BWA"
}
log.info "Exome kit:                    ${params.kit}"
if (params.panel) {
        log.info "Panel(s):                     ${params.panel}"
} else if (params.panel_intervals) {
        log.info "Panel(s):                     custom"
} else if (params.all_panels) {
        log.info "Panel(s)                      all"
}
if (params.vep) {
        log.info "Run VEP                       ${params.vep}"
}
log.info "Joint Calling                 ${params.joint_calling}"
log.info "SV Calling                    ${params.manta}"
log.info "CNV Calling		${params.cnv}"
log.info "Find expansions               ${params.expansion_hunter}"
log.info "-----------------------------------------"
log.info "Command Line:                 $workflow.commandLine"
log.info "Run name:                     ${run_name}"
if (workflow.containerEngine) {
        log.info "Container engine:             ${workflow.containerEngine}"
}
log.info "========================================="

// Read sample file
Channel.from(file(params.samples))
	.splitCsv(sep: ';', header: true)
	.map{ row-> tuple( row.IndivID,row.SampleID,row.libraryID,row.rgID,row.rgPU,row.platform,row.platform_model,row.Center,row.Date,file(row.R1),file(row.R2) ) }
	.set {  reads }

// Modules and workflows to include
include { CONVERT_BED } from "./workflows/bed/main.nf" params(params)
include { TRIM_AND_ALIGN } from "./workflows/align/main.nf" params(params)
include { VARIANT_CALLING } from "./workflows/calling/main.nf" params(params)
include { merge_gvcfs } from "./modules/deepvariant/main.nf" params(params)
include { manta } from "./modules/manta/main.nf" params(params)
include { PICARD_METRICS } from "./workflows/picard/main.nf" params(params)
include { EXPANSIONS } from "./workflows/expansions/main.nf" params(params)
include { vep } from "./modules/vep/main.nf" params(params)
include { multiqc as  multiqc_fastq ; multiqc as multiqc_library ; multiqc as multiqc_sample } from "./modules/multiqc/main.nf" params(params)
include { merge_vcf ; vcf_add_dbsnp; vcf_stats ; vcf_index } from "./modules/vcf/main.nf" params(params)
include { PANEL_QC } from "./workflows/panels/main.nf" params(params)
include { CNVKIT } from "./workflows/cnvkit/main.nf" params(params)

workflow {

	main:

		// create calling regions
		CONVERT_BED(targets)
		padded_bed = CONVERT_BED.out.bed
		bedgz = CONVERT_BED.out.bed_gz

		// align reads against genome
		TRIM_AND_ALIGN(reads)
		bam = TRIM_AND_ALIGN.out.bam
		trim_report = TRIM_AND_ALIGN.out.qc
		dedup_report = TRIM_AND_ALIGN.out.dedup_report

		// SNP + Indel calling with Deepvariant
		VARIANT_CALLING(bam,padded_bed,deepvariant_ref.collect())
		vcf = VARIANT_CALLING.out.vcf
		gvcf = VARIANT_CALLING.out.gvcf
		sample_names = VARIANT_CALLING.out.sample_names

		vcf_only = vcf.map { c,i,s,v,t -> [ v,t ] }

		// CNV Calling
		if (params.cnv) {
			CNVKIT(padded_bed,bam,cnv_cnn)
		}

		// Joint calling with GLNexus - or simple merging
		if (params.joint_calling) {
			merge_gvcfs(gvcf.collect(),padded_bed)
			merged_vcf = merge_gvcfs.out
		} else {
			merge_vcf(vcf_only.collect())
			merged_vcf = merge_vcf.out
		}
		vcf_add_dbsnp(merged_vcf)
		merged_vcf_ann = vcf_add_dbsnp.out

		// SV calling with Manta
		if (params.manta) {
			manta(bam,bedgz.collect())
			manta_vcf = manta.out[0].mix(manta.out[1],manta.out[2])
		} else {
			manta_vcf = Channel.empty()
		}

		// Expansions
		EXPANSIONS(bam,expansion_catalog)

		// QC Metrics
		PANEL_QC(bam,panels,targets)

		PICARD_METRICS(bam,targets,baits)
		bam_qc = PICARD_METRICS.out.qc_reports
		vcf_stats(vcf)
		vcf_qc = vcf_stats.out

		// Effect prediction
		vep(vcf.mix(merged_vcf_ann,manta_vcf))

		// QC
		multiqc_fastq("FastQ",trim_report.collect())
		multiqc_library("Library",dedup_report.collect())
		multiqc_sample("Sample",bam_qc.mix(vcf_qc).collect())
				
}
