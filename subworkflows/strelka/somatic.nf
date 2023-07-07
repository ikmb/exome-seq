include { STRELKA_SOMATIC } from "./../../modules/strelka/strelka_somatic"
include { BCFTOOLS_VIEW } from "./../../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../../modules/bcftools/annotate_dbsnp"

ch_versions = Channel.from([])

workflow STRELKA_SOMATIC_CALLING {

    take:
        bams
        bed
        fasta
        dbsnp

    main:

    STRELKA_SOMATIC(
        bams,
        bed,
        fasta
    )

	ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions)

    BCFTOOLS_VIEW(
        STRELKA_SOMATIC.out.vcf
    )

	ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

    BCFTOOLS_ANNOTATE_DBSNP(
        BCFTOOLS_VIEW.out.vcf.map { m,v,t ->
            [[
            patient_id: m.patient_id,
            sample_id: m.sample_id,
            variantcaller: "STRELKA"
            ],v,t]
        },
        dbsnp
        )

	ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_DBSNP.out.versions)

    emit:
	versions 	= ch_versions
    vcf 		= BCFTOOLS_ANNOTATE_DBSNP.out.vcf

}
