#!/usr/bin/env nextflow

//
// Modules and workflows to include
//

//
// Input Channels and data
//

fasta = params.genomes[ params.assembly ].fasta
fasta_fai = params.genomes[ params.assembly ].fai
fasta_gz = params.genomes[ params.assembly ].fastagz
fasta_gzfai = params.genomes[ params.assembly ].gzfai
fasta_gzi = params.genomes[ params.assembly ].gzi
dict = params.genomes[ params.assembly ].dict
dbsnp = params.genomes[ params.assembly ].dbsnp
csq_gtf = params.genomes[params.assembly].gtf
omni =  params.genomes[ params.assembly ].omni
hapmap = params.genomes[ params.assembly ].hapmap
g1k = params.genomes[ params.assembly].g1k
mills = params.genomes[ params.assembly ].mills
axiom = params.genomes[ params.assembly ].axiom

snps = [ hapmap,  omni,  dbsnp, g1k ]
indels = [ mills,  axiom ]

if (params.amplicon_bed) { ch_amplicon_bed = Channel.fromPath(file(params.amplicon_bed, checkIfExists: true)) } else { ch_amplicon_bed = Channel.from([]) }

ch_fasta = Channel.fromList( [ file(fasta , checkIfExists: true), file(fasta_fai, checkIfExits: true), file(dict, checkIfExists: true) ] )

deepvariant_ref = Channel.from( [ fasta_fai,fasta_gz,fasta_gzfai,fasta_gzi ] )

// Mapping tool and corresponding index

if (params.bwa2) { bwa_index = params.genomes[ params.assembly ].bwa2_index } else { bwa_index = fasta }

// Targets and bait file

TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

if (TARGETS==BAITS) { exit 1, "Target and bait files must not have the same name to avoid file collisions!" }

targets = Channel.from(file(TARGETS, checkIfExists: true))
baits = Channel.from(file(BAITS, checkIfExists: true))

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
sry_region  = params.sry_bed ?: params.genomes[params.assembly].sry_bed

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

// Expansion hunter references
if ('expansionhunter' in tools) {
        ecatalog = file(params.genomes[params.assembly].expansion_catalog, checkIfExists: true )
        Channel.fromPath(ecatalog)
        .ifEmpty { exit 1, "Could not find a matching ExpansionHunter catalog for this assembly" }
        .set { expansion_catalog }
} else {
        expansion_catalog = Channel.empty()
}

// CNVkit reference
if ('cnvkit' in tools) {
        cnv_ref = params.cnv_gz ?: file(params.genomes[params.assembly ].kits[params.kit].cnv_ref)
} else {
        cnv_ref = Channel.empty()
}

// Read sample file

ch_samplesheet = file(params.samples, checkIfExists: true)
ch_vcfs = Channel.from([])
ch_phased_vcfs = Channel.from([])

include { CONVERT_BED } from "./../subworkflows/bed"
include { TRIM_AND_ALIGN } from "./../subworkflows/align"
include { DV_VARIANT_CALLING } from "./../subworkflows/deepvariant"
include { GATK_VARIANT_CALLING } from "./../subworkflows/gatk"
include { STRELKA_SINGLE_CALLING } from "./../subworkflows/strelka/single"
include { STRELKA_MULTI_CALLING } from "./../subworkflows/strelka/multi"
include { GLNEXUS as MERGE_GVCFS } from "./../modules/glnexus"
include { MANTA } from "./../modules/manta.nf"
include { PICARD_METRICS } from "./../subworkflows/picard"
include { EXPANSIONS } from "./../subworkflows/expansionhunter"
include { VEP } from "./../modules/vep"
include { HAPLOSAURUS } from "./../modules/haplosaurus"
include { MULTIQC as  multiqc_fastq ; MULTIQC as multiqc_library ; MULTIQC as multiqc_sample } from "./../modules/multiqc/main"
include { BCFTOOLS_MERGE as MERGE_VCF } from "./../modules/bcftools/merge"
include { BCFTOOLS_ANNOTATE_DBSNP as VCF_ADD_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_STATS as VCF_STATS } from "./../modules/bcftools/stats"
include { TABIX as VCF_INDEX } from "./../modules/htslib/tabix"
include { PANEL_QC } from "./../subworkflows/panels"
include { CNVKIT } from "./../subworkflows/cnvkit"
include { BCFTOOLS_CSQ as CSQ } from "./../modules/bcftools/csq"
include { BCFTOOLS_CONCAT as CONCAT } from "./../modules/bcftools/concat"
include { SEX_CHECK} from "./../modules/qc/main"
include { XHLA } from "./../modules/xhla"

workflow EXOME_SEQ {

	main:

                ch_vcfs = Channel.empty()
		ch_phased_vcfs = Channel.empty()

		// create calling regions
		CONVERT_BED(targets)
		padded_bed = CONVERT_BED.out.bed
		bedgz = CONVERT_BED.out.bed_gz

		// align reads against genome
		TRIM_AND_ALIGN(
			ch_samplesheet,
			ch_amplicon_bed,
			bwa_index
		)
		bam = TRIM_AND_ALIGN.out.bam
		bam_nodedup = TRIM_AND_ALIGN.out.bam_nodedup
		trim_report = TRIM_AND_ALIGN.out.qc
		dedup_report = TRIM_AND_ALIGN.out.dedup_report
		sample_names = TRIM_AND_ALIGN.out.sample_names

		// Create a sub-set of the BAM file using a target BED file
		if ('intersect' in tools) {
			BAM_INTERSECT(
				bam,
				padded_bed
			)
		}

		// DEEPVARIANT WORKFLOW
		if ('deepvariant' in tools) {
			DV_VARIANT_CALLING(
				bam,
				padded_bed,
				deepvariant_ref.collect(),
				ch_fasta.collect(),
				dbsnp
			)
			dv_vcf = DV_VARIANT_CALLING.out.vcf
                	dv_merged_vcf = DV_VARIANT_CALLING.out.vcf_multi
			ch_vcfs = ch_vcfs.mix(dv_vcf,dv_merged_vcf)
			ch_phased_vcfs = ch_phased_vcfs.mix(
				DV_VARIANT_CALLING.out.vcf_phased_multi,
				DV_VARIANT_CALLING.out.vcf_phased_single
			)
		} else {
			dv_vcf = Channel.empty()
			dv_merged_vcf = Channel.empty()
		}

		// GATK WORKFLOW
		if ('gatk' in tools) {
			GATK_VARIANT_CALLING(
				bam,
				targets,
				TRIM_AND_ALIGN.out.metas,
				ch_fasta,
				snps,
				indels,
				dbsnp
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
					TRIM_AND_ALIGN.out.metas,
					ch_fasta
				)
				ch_vcfs = ch_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf,STRELKA_MULTI_CALLING.out.vcf_multi)
				ch_phased_vcfs = ch_phased_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf_phased_multi)
			// Call each sample individually and merge later
                        } else {
	        	        STRELKA_SINGLE_CALLING(
					bam,
					bedgz,
					sample_names,
					ch_fasta,
					dbsnp
				)
        	        	strelka_vcf = STRELKA_SINGLE_CALLING.out.vcf
				strelka_merged_vcf = STRELKA_SINGLE_CALLING.out.vcf_multi
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
		// HLA calling
		if ('xhla' in tools) {
			XHLA(bam)
		}
		// SV calling with Manta
		if ('manta' in tools) {
			MANTA(
				bam,
				bedgz.collect(),
				ch_fasta.collect()
			)
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
		PICARD_METRICS(
			bam,
			targets,
			baits,
			ch_fasta,
			dbsnp
		)
		bam_qc = PICARD_METRICS.out.qc_reports
		VCF_STATS(ch_vcfs)
		vcf_qc = VCF_STATS.out.stats

		// Coverage of SRY gene
		SEX_CHECK(
			bam.map { m,b,i -> tuple(b,i) }.collect(),
			sry_region
		)

		// Effect prediction
		if ('vep' in tools) {
			VEP(ch_vcfs,ch_fasta)
		}
		if ('csq' in tools) {
			CSQ(ch_phased_vcfs,ch_fasta)
		}
		if ('haplosaurus' in tools) {
			HAPLOSAURUS(ch_phased_vcfs,ch_fasta)
		}

		// QC Reports
		multiqc_fastq("FastQ",trim_report.collect())
		multiqc_library("Library",dedup_report.collect())
		multiqc_sample("Sample",bam_qc.mix(vcf_qc,SEX_CHECK.out.yaml).collect())


	emit:
		qc = multiqc_sample.out.report.toList()

}
