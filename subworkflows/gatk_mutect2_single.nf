include { GATK_MUTECT2 } from "./../modules/gatk/mutect2"
include { GATK_FILTER_MUTECT_CALLS } from "./../modules/gatk/filter_mutect_calls"
include { BCFTOOLS_VIEW } from "./../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_ANNOTATE } from "./../modules/bcftools/annotate"
include { GATK_SPLITINTERVALS } from './../modules/gatk/splitintervals'
include { GATK_LEARN_READ_ORIENTATION_MODEL } from "./../modules/gatk/learn_read_orientation_model"
include { GATK_GET_PILEUP_SUMMARIES } from "./../modules/gatk/get_pileup_summaries"
include { GATK_CALCULATE_CONTAMINATION } from "./../modules/gatk/calculate_contamination"
include { PICARD_MERGEVCFS } from "./../modules/picard/mergevcfs"

ch_versions = Channel.from([])
ch_vcfs 	= Channel.from([])

workflow GATK_MUTECT2_SINGLE {

    take:
    bam
    targets
    fasta
    dbsnp
    mutect_normals
    mutect_normals_tbi
	
    main:

 // Make split targets
    GATK_SPLITINTERVALS(
        targets
    )
    GATK_SPLITINTERVALS.out.intervals.flatMap { i ->
        i.collect { file(it) }
    }.set { targets_split }

    GATK_MUTECT2(
        bam,
        targets_split,
        fasta,
        mutect_normals,
        mutect_normals_tbi
    )

    ch_versions = ch_versions.mix(GATK_MUTECT2.out.versions)

    PICARD_MERGEVCFS(
        GATK_MUTECT2.out.vcf.groupTuple()
    )

    GATK_LEARN_READ_ORIENTATION_MODEL(
        GATK_MUTECT2.out.f1r2.groupTuple()
    )

    ch_versions = ch_versions.mix(GATK_LEARN_READ_ORIENTATION_MODEL.out.versions)

    GATK_GET_PILEUP_SUMMARIES(
        bam,
        targets,
        fasta
    )

    ch_versions = ch_versions.mix(GATK_GET_PILEUP_SUMMARIES.out.versions)

    GATK_CALCULATE_CONTAMINATION(
        GATK_GET_PILEUP_SUMMARIES.out.table
    )

    ch_versions = ch_versions.mix(GATK_CALCULATE_CONTAMINATION.out.versions)

    ch_mutect = PICARD_MERGEVCFS.out.vcf.join(
        GATK_LEARN_READ_ORIENTATION_MODEL.out.model
    ).join(GATK_CALCULATE_CONTAMINATION.out.table)

    GATK_FILTER_MUTECT_CALLS(
        ch_mutect,
        fasta
    )

    ch_versions = ch_versions.mix(GATK_FILTER_MUTECT_CALLS.out.versions)

    ch_vcfs = ch_vcfs.mix(GATK_FILTER_MUTECT_CALLS.out.vcf)

    BCFTOOLS_VIEW(ch_vcfs)

    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    BCFTOOLS_ANNOTATE_DBSNP(
        BCFTOOLS_VIEW.out.vcf.map { meta,v,t ->
            def s_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "MUTECT2" ]
            tuple(s_meta,v,t)
        },
        dbsnp
    )
		
    ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_DBSNP.out.versions)

    BCFTOOLS_ANNOTATE(
        BCFTOOLS_ANNOTATE_DBSNP.out.vcf
    )

    emit:
    versions 	= ch_versions
    vcf 	= BCFTOOLS_ANNOTATE.out.vcf
		
}

