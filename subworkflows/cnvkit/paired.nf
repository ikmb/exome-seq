include { GUNZIP } from "./../../modules/gunzip"
include { CNVKIT_BATCH_PAIRED } from "./../../modules/cnvkit/batch_paired"
include { CNVKIT_EXPORT } from "./../../modules/cnvkit/export"

ch_versions = Channel.from([])

workflow CNVKIT_PAIRED{

    take:
        bam
        targets
        fasta

    main:

        CNVKIT_BATCH_PAIRED(
            bam,
            targets.collect(),
            fasta.collect()
        )
    
        ch_versions = ch_versions.mix(CNVKIT_BATCH_PAIRED.out.versions)

        CNVKIT_EXPORT(
            CNVKIT_BATCH_PAIRED.out.cns
        )
        
    emit:
    versions = ch_versions
}
