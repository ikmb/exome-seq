include { PANEL_COVERAGE } from "./../modules/picard/panel_coverage"
include { MULTIQC_PANEL } from "./../modules/multiqc/main.nf"

ch_versions = Channel.from([])

workflow PANEL_QC {

    take:
        bam
        panels
        targets

    main:
        PANEL_COVERAGE(
            bam.combine(panels),
            targets
        )

        ch_versions = ch_versions.mix(PANEL_COVERAGE.out.versions)

        MULTIQC_PANEL(
            PANEL_COVERAGE.out.coverage.groupTuple()
        )
    

    emit:
        versions = ch_versions
        qc = MULTIQC_PANEL.out

}
    
