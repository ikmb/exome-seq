include { GUNZIP } from "./../../modules/gunzip"
include { CNVKIT_BATCH } from "./../../modules/cnvkit/batch"
include { CNVKIT_EXPORT } from "./../../modules/cnvkit/export"

ch_versions = Channel.from([])

workflow CNVKIT_SINGLE {

    take:
        bam
        cnn_gz

    main:

        CNVKIT_BATCH(
            bam,
            cnn_gz.collect()
        )
    
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

        CNVKIT_EXPORT(
            CNVKIT_BATCH.out.cns
        )
        
    emit:
    versions = ch_versions
}
