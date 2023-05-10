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

        ch_cns = Channel.from([])

        // Results can be one or multiple CNS files, depending on how many tumor samples were used per patient. Need to normalize into meta,cns emission
        ch_cns = ch_cns.mix( CNVKIT_BATCH_PAIRED.out.cns.filter { m,cns -> !(cns instanceof List) } )
        ch_cns = ch_cns.mix(
            CNVKIT_BATCH_PAIRED.out.cns.filter { m,cns -> cns instanceof List }.flatMap { m,cns ->
                cns.collect{ [ m,file(it) ] }
            }
        )

        // this will emit multiple call.cns files so we need to split this into singletons again
        CNVKIT_EXPORT(
            ch_cns
        )
        

    emit:
    versions = ch_versions
}
