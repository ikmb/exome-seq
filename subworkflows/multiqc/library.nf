include { MULTIQC } from './../../modules/multiqc/main.nf'

workflow multiqc_library {
    take:
    reports


    main:
    MULTIQC("Library", reports)

    emit:
    report = MULTIQC.out.report
}