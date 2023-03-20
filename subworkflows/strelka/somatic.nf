include { STRELKA_SOMATIC } from "./../../modules/strelka/strelka_somatic"
include { BCFTOOLS_VIEW } from "./../../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../../modules/bcftools/annotate_dbsnp"

workflow STRELKA_SOMATIC_CALLING {

    take:
        bams
        bed
        fasta
        dbsnp

    main:

        STRELKA_SOMATIC(
            bams,
            bed.collect(),
            fasta.collect()
        )

        BCFTOOLS_VIEW(
            STRELKA_SOMATIC.out.vcf
        )

        BCFTOOLS_ANNOTATE_DBSNP(
            BCFTOOLS_VIEW.out.vcf.map { m,v,t ->
                [[
                    patient_id: m.patient_id,
                    sample_id: m.sample_id,
                    variantcaller: "STRELKA"
                ],v,t]
            },
            dbsnp.collect()
        )

    emit:
        vcf = BCFTOOLS_ANNOTATE_DBSNP.out.vcf

}
