include { GUNZIP } from "./../../modules/gunzip"
include { CNVKIT_BATCH } from "./../../modules/cnvkit/batch"
include { CNVKIT_EXPORT } from "./../../modules/cnvkit/export"

ch_versions = Channel.from([])

workflow CNVKIT_SINGLE {

    take:
        bam
        cnn_gz
        fasta

    main:
        
        GUNZIP(
            Channel.fromPath(cnn_gz).map { f -> 
                tuple([ id: f.getBaseName() ], f)
            }
        )

        CNVKIT_BATCH(
            bam,
            GUNZIP.out.decompressed.collect(),
            fasta.collect()
        )
    
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

        CNVKIT_EXPORT(
            CNVKIT_BATCH.out.cns
        )
        
    emit:
    versions = ch_versions
}
