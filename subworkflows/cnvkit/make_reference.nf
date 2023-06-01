include { CNVKIT_MAKE_FLAT_REF } from "../../modules/cnvkit/make_flat_ref"
include { GZIP } from "../../modules/gzip"

workflow CNVKIT_MAKE_REFERENCE {

    take:
    bams
    bed
    fasta

    main:
 
    CNVKIT_MAKE_FLAT_REF(
        bed.collect(),
        fasta.collect()
    )

    cnv_ref = CNVKIT_MAKE_FLAT_REF.out.cnn


    emit:
    cnn = cnv_ref

}

