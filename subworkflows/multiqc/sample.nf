include { MULTIQC } from './../../modules/multiqc/main.nf'

workflow multiqc_sample {
    take:
    reports

    main:
    MULTIQC("Sample", reports)

    emit:
    report = MULTIQC.out.report
}