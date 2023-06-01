#!/usr/bin/env nextflow

//
// Modules and workflows to include
//

// ************************************
// Input Channels and data
// ************************************

ch_dbsnp        = Channel.fromPath(params.dbsnp)
ch_dbsnp_tbi    = Channel.fromPath(params.dbsnp + ".tbi")
ch_hapmap       = Channel.fromPath(params.hapmap)
ch_hapmap_tbi   = Channel.fromPath(params.hapmap + ".tbi")
ch_omni         = Channel.fromPath(params.omni)
ch_omni_tbi     = Channel.fromPath(params.omni + ".tbi")
ch_mills        = Channel.fromPath(params.mills)
ch_mills_tbi    = Channel.fromPath(params.mills + ".tbi")
ch_g1k          = Channel.fromPath(params.g1k)
ch_g1k_tbi      = Channel.fromPath(params.g1k + ".tbi")
ch_axiom        = Channel.fromPath(params.axiom)
ch_axiom_tbi    = Channel.fromPath(params.axiom + ".tbi")

// ************************************
// combine all SNPs, for GATK calibration
// ************************************

ch_known_snps		= ch_dbsnp.mix(ch_hapmap, ch_omni, ch_g1k).collect()
ch_known_snps_tbi 	= ch_dbsnp_tbi.mix(ch_hapmap_tbi, ch_omni_tbi, ch_g1k_tbi).collect()

// ************************************
// combine all INDELs, for GATK calibration
// ************************************
ch_known_indels     = ch_mills.mix(ch_axiom).collect()
ch_known_indels_tbi = ch_mills_tbi.mix(ch_axiom_tbi).collect()

ch_dbsnp_combined   = Channel.fromList( [ params.dbsnp, params.dbsnp + ".tbi" ] ).collect()

// ************************************
// Provide a BED file with amplicon locations
// ************************************
if (params.amplicon_bed) { ch_amplicon_bed = Channel.fromPath(file(params.amplicon_bed, checkIfExists: true)) } else { ch_amplicon_bed = Channel.from([]) }

// *************************************
// A GTF file for protein-level effect prediction
// *************************************
ch_gtf              = Channel.fromPath(params.csq_gtf)

// *************************************
// The reference genome with relevant helper files
// *************************************
ch_fasta = Channel.fromList( [ file(params.fasta , checkIfExists: true), file(params.fasta_fai, checkIfExits: true), file(params.dict, checkIfExists: true) ] ).collect()

// ************************************
// Mapping tool and corresponding index
// ************************************
if (!params.aligners_allowed.contains(params.aligner)) {
    exit 1, "Invalid aligner specified - valid options are: ${params.aligners_allowed.join(',')}"
}

// ************************************
// Choose your aligner
// ************************************
if (params.aligner == "dragmap") {
    genome_index        = params.dragmap_index
} else if (params.aligner == "bwa2") {
        genome_index    = params.bwa2_index
} else {
    genome_index        = params.bwa_index
}

// ************************************
// CNVkit reference
// ************************************
cnv_ref = false
if (params.skip_cnv_gz) {
    ch_cnv_gz 		= Channel.value([])
} else if (params.cnv_gz) {
    cnv_ref             = true
    ch_cnv_gz 		= Channel.fromPath(params.cnv_gz)
} else if ( params.kit && params.genomes[ params.assembly ].kits[ params.kit ].cnv_ref) { 
    cnv_ref             = true
    ch_cnv_gz 		= Channel.fromPath(params.genomes[ params.assembly ].kits[ params.kit ].cnv_ref)
} else { 
    ch_cnv_gz 		= Channel.value([])
}

// ************************************
// Targets and bait file
// ************************************

if (!params.kit && !params.targets || !params.kit && !params.baits) {
    exit 1, "No exome kit or custom targets (--targets) and/or baits (--baits) file specified!"
}
TARGETS = params.targets ?: params.genomes[params.assembly].kits[ params.kit ].targets
BAITS = params.baits ?: params.genomes[params.assembly].kits[ params.kit ].baits

if (TARGETS==BAITS) { exit 1, "Target and bait files must not have the same name to avoid file collisions!" }

if (!TARGETS || !BAITS) { exit 1, "No kit or user-supplied calling intervals found - cannot proceed!" }

targets = Channel.from(file(TARGETS, checkIfExists: true)).collect()
baits = Channel.from(file(BAITS, checkIfExists: true)).collect()

// ************************************
//PANEL COVERAGE - pick the correct panel for reporting
// ************************************

if (params.panel) {
    panel           = params.genomes[params.assembly].panels[params.panel].intervals
    panels          = Channel.fromPath(panel)
} else if (params.panel_intervals) {
    Channel.fromPath(params.panel_intervals)
    .ifEmpty { exit 1; "Could not find the specified gene panel (--panel_intervals)" }
    .set { panels }
} else if (params.all_panels) {
    panel_list 		= []
    panel_names 	= params.genomes[params.assembly].panels.keySet()
    panel_names.each {
        interval 	= params.genomes[params.assembly].panels[it].intervals
        panel_list << file(interval)
    }
    panels = Channel.fromList(panel_list)
} else {
    panels = Channel.empty()
}

// ************************************
// A single exon only covered in male samples - simple sex check
// ************************************
sry_region  = params.sry_bed ?: params.genomes[params.assembly].sry_bed

// ************************************
// List of tools to run
// ************************************
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []

// ************************************
// Expansion hunter references
// ************************************
if ('expansionhunter' in tools) {
        ecatalog = file(params.genomes[params.assembly].expansion_catalog, checkIfExists: true )
        Channel.fromPath(ecatalog)
        .ifEmpty { exit 1, "Could not find a matching ExpansionHunter catalog for this assembly" }
        .set { expansion_catalog }
} else {
        expansion_catalog = Channel.empty()
}

// **********************************
// Mutect2 panel of normals, if any
// **********************************
if ('mutect2' in tools) {

    // Skip PONs no matter what
    if (params.skip_mutect_pon) {
        ch_mutect_pon 		= Channel.value([])
        ch_mutect_pon_tbi 	= Channel.value([])
    // Use a user-supplied PON
    }else if (params.mutect_normals) {
        ch_mutect_pon 		= Channel.fromPath(file(params.mutect_normals), checkIfExists: true).collect()
        ch_mutect_pon_tbi 	= Channel.fromPath(file(params.mutect_normals + ".tbi"), checkIfExists: true).collect()
    // See if there is a pre-configured PON
    } else if ( params.kit && params.genomes[params.assembly].kits[params.kit].mutect_pon) {
        ch_mutect_pon 		= Channel.fromPath(file(params.genomes[params.assembly].kits[params.kit].mutect_pon), checkIfExists: true).collect()
        ch_mutect_pon_tbi 	= Channel.fromPath(file(params.genomes[params.assembly].kits[params.kit].mutect_pon + ".tbi"), checkIfExists: true).collect()
    // Or else return empty channels
    } else {
        ch_mutect_pon 		= Channel.value([])
        ch_mutect_pon_tbi 	= Channel.value([])
    }
} else {
    ch_mutect_pon           = Channel.value([])
    ch_mutect_pon_tbi 		= Channel.value([])
}

// ************************************
// Read sample file
// ************************************

ch_samplesheet 			= Channel.fromPath(params.samples)

// ************************************
// set optional channels
// ************************************

ch_vcfs                 = Channel.from([])
ch_phased_vcfs 			= Channel.from([])
ch_recal_bam 			= Channel.from([])
ch_manta_indels 		= Channel.from([])
ch_multiqc_files 		= Channel.from([])
ch_versions 			= Channel.from([])
ch_manta_vcfs 			= Channel.from([])
ch_manta_indels_paired  = Channel.from([])
ch_bam_normal			= Channel.from([])

// ************************************
// import subworkflows and modules
// ************************************

include { CONVERT_BED } from "./../subworkflows/bed"
include { TRIM_AND_ALIGN } from "./../subworkflows/align"
include { DV_VARIANT_CALLING } from "./../subworkflows/deepvariant"
include { GATK_VARIANT_CALLING } from "./../subworkflows/gatk_variant_calling"
include { GATK_BAM_RECAL } from "./../subworkflows/gatk_bqsr"
include { GATK_MUTECT2_SINGLE } from "./../subworkflows/gatk_mutect2_single"
include { GATK_MUTECT2_PAIRED } from "./../subworkflows/gatk_mutect2_paired"
include { STRELKA_SINGLE_CALLING } from "./../subworkflows/strelka/single"
include { STRELKA_MULTI_CALLING } from "./../subworkflows/strelka/multi"
include { STRELKA_SOMATIC_CALLING } from "./../subworkflows/strelka/somatic"
include { GLNEXUS as MERGE_GVCFS } from "./../modules/glnexus"
include { MANTA_NORMAL } from "./../modules/manta/normal"
include { MANTA_TUMOR } from "./../modules/manta/tumor"
include { MANTA_PAIRED } from "./../modules/manta/paired"
include { PICARD_METRICS } from "./../subworkflows/picard"
include { EXPANSIONS } from "./../subworkflows/expansionhunter"
include { VEP_VEP as VEP } from "./../modules/vep/vep"
include { HAPLOSAURUS } from "./../modules/haplosaurus"
include { MULTIQC as  multiqc_fastq ; MULTIQC as multiqc_library ; MULTIQC as multiqc_sample } from "./../modules/multiqc/main"
include { BCFTOOLS_MERGE as MERGE_VCF } from "./../modules/bcftools/merge"
include { BCFTOOLS_ANNOTATE_DBSNP as VCF_ADD_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_STATS as VCF_STATS } from "./../modules/bcftools/stats"
include { TABIX as VCF_INDEX } from "./../modules/htslib/tabix"
include { PANEL_QC } from "./../subworkflows/panels"
include { BCFTOOLS_CSQ as CSQ } from "./../modules/bcftools/csq"
include { BCFTOOLS_CONCAT as CONCAT } from "./../modules/bcftools/concat"
include { SEX_CHECK} from "./../modules/qc/main"
include { XHLA } from "./../modules/xhla"
include { CNVKIT_SINGLE } from "./../subworkflows/cnvkit/single"
include { CNVKIT_PAIRED } from "./../subworkflows/cnvkit/paired"
include { CNVKIT_MAKE_REFERENCE } from "./../subworkflows/cnvkit/make_reference"
include { VALIDATE_SAMPLESHEET } from "./../modules/validate_samplesheet"
include { CUSTOM_DUMPSOFTWAREVERSIONS } from "./../modules/custom/dumpsoftwareversions/main"

// Start the main workflow
workflow EXOME_SEQ {

    main:

    // create calling regions
    CONVERT_BED(
        targets,
        ch_fasta
    )
    padded_bed  = CONVERT_BED.out.bed_padded.collect()
    bedgz		= CONVERT_BED.out.bed_gz.collect()
    bed 		= CONVERT_BED.out.bed.collect()

    // Make sure the format of the samplesheet is correct
    VALIDATE_SAMPLESHEET(
        ch_samplesheet
    )

    // Trim and align reads against genome
    TRIM_AND_ALIGN(
        VALIDATE_SAMPLESHEET.out.csv,
        ch_amplicon_bed,
        genome_index,
        ch_fasta
    )
    ch_bam		= TRIM_AND_ALIGN.out.bam
    ch_bam_nodedup	= TRIM_AND_ALIGN.out.bam_nodedup
    trim_report		= TRIM_AND_ALIGN.out.qc
    dedup_report	= TRIM_AND_ALIGN.out.dedup_report
    sample_names	= TRIM_AND_ALIGN.out.sample_names

    ch_versions		= ch_versions.mix(TRIM_AND_ALIGN.out.versions)

    // Create a sub-set of the BAM file using a target BED file
    if ('intersect' in tools) {
        BAM_INTERSECT(
            ch_bam,
            padded_bed
        )
    }

    // DEEPVARIANT WORKFLOW
    if ('deepvariant' in tools) {
        DV_VARIANT_CALLING(
            ch_bam,
            padded_bed,
            ch_fasta,
            ch_dbsnp_combined
        )
        dv_vcf            = DV_VARIANT_CALLING.out.vcf
        dv_merged_vcf     = DV_VARIANT_CALLING.out.vcf_multi
        ch_vcfs           = ch_vcfs.mix(dv_vcf,dv_merged_vcf)

        ch_phased_vcfs = ch_phased_vcfs.mix(
            DV_VARIANT_CALLING.out.vcf_phased_multi,
            DV_VARIANT_CALLING.out.vcf_phased_single
        )
            
        ch_versions = ch_versions.mix(DV_VARIANT_CALLING.out.versions)
    } else {
        dv_vcf = Channel.empty()
        dv_merged_vcf = Channel.empty()
    }

    // Make GATK-compliant BAM file
    if ('gatk' in tools || 'mutect2' in tools) {
        GATK_BAM_RECAL(
            ch_bam,
            targets,
            ch_fasta,
            ch_known_snps,
            ch_known_snps_tbi,
            ch_known_indels,
            ch_known_indels_tbi
        )
        ch_recal_bam    = ch_recal_bam.mix(GATK_BAM_RECAL.out.bam)
        ch_versions     = ch_versions.mix(GATK_BAM_RECAL.out.versions)
    }
                
    // GATK HAPLOTYPECALLER WORKFLOW
    if ('gatk' in tools) {
        GATK_VARIANT_CALLING(
            ch_recal_bam,
            targets,
            ch_fasta,
            ch_known_snps,
            ch_known_snps_tbi,
            ch_known_indels,
            ch_known_indels_tbi    
        )
        gatk_vcf        = GATK_VARIANT_CALLING.out.vcf
        gatk_merged_vcf = GATK_VARIANT_CALLING.out.vcf_multi
        ch_vcfs         = ch_vcfs.mix(GATK_VARIANT_CALLING.out.vcf_multi)
        ch_versions     = ch_versions.mix(GATK_VARIANT_CALLING.out.versions)
            
    } else {
        gatk_vcf = Channel.empty()
        gatk_merged_vcf = Channel.empty()
    }

    // MUTECT2 workflow - only works with tumor or tumor-normal pairs. 
    if ('mutect2' in tools) {

        ch_recal_bam.branch { m,b,i ->
            normal: m.status == 0
            tumor: m.status == 1
        }.set { ch_recal_bam_status }
            
        // Fetch tumor and normal samples and group by patient ID into channel for paired calling (if any)
        ch_recal_bam_normal                         = ch_recal_bam_status.normal
        ch_recal_bam_tumor                          = ch_recal_bam_status.tumor

        ch_recal_bam_normal_cross                   = ch_recal_bam_normal.map { m,b,i -> [ m.patient_id,m,b,i] }
        ch_recal_bam_tumor_cross                    = ch_recal_bam_tumor.map { m,b,i -> [ m.patient_id,m,b,i] }

        // all tumor samples belonging to the same patient
        ch_recal_bam_tumor_grouped                  = ch_recal_bam_tumor_cross.groupTuple()
        ch_recal_bam_tumor_grouped_joined           = ch_recal_bam_tumor_grouped.join(ch_recal_bam_normal, remainder: true)
        ch_recal_bam_tumor_grouped_joined_filtered  = ch_recal_bam_tumor_grouped_joined.filter{ it -> !(it.last()) }
        ch_recal_bam_tumor_grouped_tumor_only       = ch_recal_bam_tumor_grouped_joined_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }
        
        // combining each normal with all matched tumor samples for joint analysis
        ch_recal_bam_normal_cross_joined            = ch_recal_bam_normal_cross.join(ch_recal_bam_tumor_grouped)
        ch_recal_bam_normal_cross_joined_filtered   = ch_recal_bam_normal_cross_joined.filter{ it -> it.last() }

        ch_recal_bam_normal_cross_joined_filtered.map { i,m,nb,nbi,mt,tb,tbi ->
            [[
            patient_id: m.patient_id,
            normal_id: m.sample_id,
            tumor_id: "ALL_TUMOR_SAMPLES",
            sample_id: "${m.sample_id}_vs_all_tumors".toString()
            ],nb,nbi,tb,tbi]
        }.set { ch_recal_bam_normal_grouped_tumor }

	ch_recal_bam_normal_grouped_tumor.view()

        // combining each normal sample with each tumor sample for pair-wise analysis
        ch_recal_bam_tumor_joined 			= ch_recal_bam_tumor_cross.join(ch_recal_bam_normal_cross, remainder: true)
        ch_recal_bam_tumor_joined_filtered              = ch_recal_bam_tumor_joined.filter{ it ->  !(it.last()) }
        ch_recal_bam_tumor_only 			= ch_recal_bam_tumor_joined_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

        ch_recal_bam_normal_cross.cross(ch_recal_bam_tumor_cross).map { normal,tumor ->
            [[
                patient_id: normal[0],
                normal_id: "${normal[1].sample_id}",
                tumor_id: "${tumor[1].sample_id}",
                sample_id: "${tumor[1].sample_id}_vs_${normal[1].sample_id}".toString()    
            ],normal[2],normal[3],tumor[2],tumor[3]]
        
        }.set { ch_recal_bam_calling_pair }

        // Variant calling for paired tumor-normal samples
        GATK_MUTECT2_PAIRED(
            ch_recal_bam_normal_grouped_tumor,
            targets,
            ch_fasta,
            ch_dbsnp_combined,
            ch_mutect_pon,
            ch_mutect_pon_tbi
        )

        ch_vcfs         = ch_vcfs.mix(GATK_MUTECT2_PAIRED.out.vcf)
        ch_versions     = ch_versions.mix(GATK_MUTECT2_PAIRED.out.versions)

        // Variant calling for tumor-only samples
        GATK_MUTECT2_SINGLE(
            ch_recal_bam_tumor_only,
            targets,
            ch_fasta,
            ch_dbsnp_combined,
            ch_mutect_pon,
            ch_mutect_pon_tbi
        )

        ch_vcfs         = ch_vcfs.mix(GATK_MUTECT2_SINGLE.out.vcf)
        ch_versions     = ch_versions.mix(GATK_MUTECT2_SINGLE.out.versions)
    }

    // *********************************
    // Create channels for different calling strategies re: tumor/normal
    // *********************************

    ch_bam.branch { m,b,i ->
        normal: m.status == 0
        tumor: m.status == 1
    }.set { ch_bam_status }

    // Fetch tumor and normal samples and group by patient ID into channel for paired calling (if any)
    ch_bam_normal           = ch_bam_normal.mix(ch_bam_status.normal)
    ch_bam_tumor            = ch_bam_status.tumor

    ch_bam_normal_cross     = ch_bam_normal.map { m,b,i -> [ m.patient_id,m,b,i] }
    ch_bam_tumor_cross      = ch_bam_tumor.map { m,b,i -> [ m.patient_id,m,b,i] }

    // Find and group unmatched tumor samples
    ch_bam_tumor_cross_grouped                      = ch_bam_tumor_cross.groupTuple()
    ch_bam_tumor_cross_grouped_joined               = ch_bam_tumor_cross_grouped.join(ch_bam_normal_cross, remainder: true)
    ch_bam_tumor_cross_grouped_joined_filtered      = ch_bam_tumor_cross_grouped_joined.filter{ it -> !(it.last()) }
    ch_bam_tumor_cross_grouped_tumor_only           = ch_bam_tumor_cross_grouped_joined_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

    // combine each normal sample with all matched tumor samples for joint analysis
    ch_bam_normal_cross_joined                      = ch_bam_normal_cross.join(ch_bam_tumor_cross_grouped)
    ch_bam_normal_cross_joined_filtered             = ch_bam_normal_cross_joined.filter{ it -> it.last() }

    ch_bam_normal_cross_joined_filtered.map { i,m,nb,nbi,mt,tb,tbi ->
        [[
        patient_id: m.patient_id,
        normal_id: m.sample_id,
        tumor_id: "ALL_TUMOR_SAMPLES",
        sample_id: "${m.sample_id}_vs_all_tumors".toString()
        ],nb,nbi,tb,tbi]
    }.set { ch_bam_normal_grouped_tumor }
    
    // combine each normal sample with each matched tumor sample for _pairwise_ analysis
    ch_bam_normal_cross.cross(ch_bam_tumor_cross).map { normal,tumor ->
        [[
        patient_id: normal[0],
        normal_id: "${normal[1].sample_id}",
        tumor_id: "${tumor[1].sample_id}",
        sample_id: "${tumor[1].sample_id}_vs_${normal[1].sample_id}".toString()
        ],normal[2],normal[3],tumor[2],tumor[3]]
    }.set { ch_bam_calling_pair }

    // *********************
    // SV calling with Manta
    // *********************

    if ('manta' in tools || 'strelka' in tools ) {
        
        MANTA_NORMAL(
            ch_bam_normal,
            bedgz.collect(),
            ch_fasta.collect()
        )

        ch_manta_indels = ch_manta_indels.mix(MANTA_NORMAL.out.small_indels)
        ch_manta_vcfs = ch_manta_vcfs.mix(MANTA_NORMAL.out.diploid_sv)

        ch_versions = ch_versions.mix(MANTA_NORMAL.out.versions)

        MANTA_PAIRED(
            ch_bam_calling_pair,
            bedgz.collect(),
            ch_fasta.collect()
        )

        ch_manta_indels_paired = ch_manta_indels_paired.mix(MANTA_PAIRED.out.small_indels)

        ch_versions = ch_versions.mix(MANTA_PAIRED.out.versions)
        ch_manta_vcfs = ch_manta_vcfs.mix(MANTA_PAIRED.out.diploid_sv,MANTA_PAIRED.out.somatic_sv)

        MANTA_TUMOR(
            ch_bam_tumor_cross_grouped_tumor_only,
            bedgz.collect(),
            ch_fasta.collect()
        )
        ch_manta_indels = ch_manta_indels.mix(MANTA_TUMOR.out.small_indels)

        ch_versions = ch_versions.mix(MANTA_TUMOR.out.versions)
        ch_manta_vcfs = ch_manta_vcfs.mix(MANTA_TUMOR.out.tumor_sv)

    } else {
        ch_manta_vcfs = Channel.empty()
    }

    // STRELKA WORKFLOW
    if ('strelka' in tools) {

        // Join Manta normal Indel calls to tumor-normal pair via the normal_id/sample_id
        ch_bam_calling_pair_manta = ch_bam_calling_pair.join(ch_manta_indels_paired)

        // Tumor-normal pairs only with Manta
        STRELKA_SOMATIC_CALLING(
            ch_bam_calling_pair_manta,
            bedgz,
            ch_fasta,
            ch_dbsnp_combined
        )
        ch_vcfs = ch_vcfs.mix(STRELKA_SOMATIC_CALLING.out.vcf)
        ch_versions = ch_versions.mix(STRELKA_SOMATIC_CALLING.out.versions)

        // Call all samples together
        if (params.joint_calling) {
            STRELKA_MULTI_CALLING(
                ch_bam.map {m,b,i -> [ b,i ] },
                bedgz,
                TRIM_AND_ALIGN.out.metas,
                ch_fasta
            )
            ch_vcfs = ch_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf,STRELKA_MULTI_CALLING.out.vcf_multi)
            ch_phased_vcfs = ch_phased_vcfs.mix(STRELKA_MULTI_CALLING.out.vcf_phased_multi)
            ch_versions = ch_versions.mix(STRELKA_MULTI_CALLING.out.versions)

        // Call each sample individually and merge later
        } else {

            ch_bam_manta_single = ch_bam.join(ch_manta_indels )

            STRELKA_SINGLE_CALLING(
                ch_bam_manta_single,
                bedgz,
                sample_names,
                ch_fasta,
                ch_dbsnp_combined
            )
            strelka_vcf         = STRELKA_SINGLE_CALLING.out.vcf
            strelka_merged_vcf  = STRELKA_SINGLE_CALLING.out.vcf_multi
            ch_vcfs             = ch_vcfs.mix(strelka_vcf,strelka_merged_vcf)
            ch_versions         = ch_versions.mix(STRELKA_SINGLE_CALLING.out.versions)
        }
    } else {
        strelka_vcf = Channel.empty()
        strelka_merged_vcf = Channel.empty()
    }
        
    // HLA calling
    if ('xhla' in tools) {
        XHLA(ch_bam)
    }

    // CNV calling
    if ('cnvkit' in tools) {

        if (cnv_ref) {
            cnv_cnn_ref = ch_cnv_gz
        } else {            
            CNVKIT_MAKE_REFERENCE(
                ch_bam_normal,
                bed,
                ch_fasta
            )
            cnv_cnn_ref = CNVKIT_MAKE_REFERENCE.out.cnn
        }

        CNVKIT_SINGLE(
            ch_bam,
            cnv_cnn_ref
        )

        CNVKIT_PAIRED(
            ch_bam_normal_grouped_tumor,
            bed,
            ch_fasta
        )

        ch_versions = ch_versions.mix(CNVKIT_PAIRED.out.versions)
    }
        
    // Expansions
    if ('expansionhunter' in tools) {
        EXPANSIONS(
            ch_bam,
            expansion_catalog.collect()
        )
        ch_versions = ch_versions.mix(EXPANSIONS.out.versions)
    }

    // QC Metrics
    PANEL_QC(
        ch_bam,
        panels,
        targets
    )
    PICARD_METRICS(
        ch_bam,
        targets,
        baits,
        ch_fasta,
        ch_dbsnp_combined
    )
    bam_qc = PICARD_METRICS.out.qc_reports
    VCF_STATS(ch_vcfs)
    vcf_qc = VCF_STATS.out.stats

    ch_versions = ch_versions.mix(PICARD_METRICS.out.versions)

    // Coverage of SRY gene
    SEX_CHECK(
        ch_bam.map { m,b,i -> tuple(b,i) }.collect(),
        sry_region
    )

    // Effect prediction
    if ('vep' in tools) {
        VEP(
            ch_vcfs,
            ch_fasta.collect()
        )
            //ch_versions = ch_versions.mix(VEP.out.versions)
    }
    if ('csq' in tools) {
        CSQ(
            ch_phased_vcfs,
            ch_fasta.collect(),
            ch_gtf.collect()
        )

        //ch_versions = ch_versions.mix(CSQ.out.versions)

    }
    if ('haplosaurus' in tools) {
        HAPLOSAURUS(
            ch_phased_vcfs,
            ch_fasta.collect()
        )

        ch_versions = ch_versions.mix(HAPLOSAURUS.out.versions)

    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // QC Reports
    multiqc_fastq(
        "FastQ",
        trim_report.collect()
    )
    multiqc_library(
        "Library",
        dedup_report.collect()
    )
    multiqc_sample(
        "Sample",
        bam_qc.mix(
            vcf_qc,
            SEX_CHECK.out.yaml,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml
        ).collect()
    )

    emit:
    qc = multiqc_sample.out.report

}
