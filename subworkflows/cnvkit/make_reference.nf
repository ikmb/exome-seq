include { CNVKIT_MAKE_REF } from "../../modules/cnvkit/make_ref" 
include { GZIP } from "../../modules/gzip"

workflow CNVKIT_MAKE_REFERENCE {

    take:
    bams
    bed
    fasta

    main:

    CNVKIT_MAKE_REF(
        bams.map { m,b,i -> 
           [ b,i ]
        }.collect(),
        bed.collect(),
        fasta.collect(),
        bams.count()
    )

    emit:
    cnn = CNVKIT_MAKE_REF.out.cnn

}

