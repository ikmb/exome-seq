include { PANEL_COVERAGE } from "./../modules/picard/panel_coverage"
include { MULTIQC_PANEL } from "./../modules/multiqc/main.nf"

workflow PANEL_QC {

	take:
		bam
		panels
		targets

	main:
		PANEL_COVERAGE(
			bam.combine(panels),
			targets.collect()
		)
		MULTIQC_PANEL(
			PANEL_COVERAGE.out.coverage.groupTuple()
		)
	

	emit:
		qc = MULTIQC_PANEL.out

}
	
