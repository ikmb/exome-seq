include { MULTIQC } from './../../modules/multiqc/main.nf'

workflow multiqc_fastq {
    take:
    reports

    main:
    MULTIQC("FastQ", reports)

    emit:
    report = MULTIQC.out.report
}