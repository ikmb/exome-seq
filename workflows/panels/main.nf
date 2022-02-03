include {  panel_coverage } from "./../../modules/picard/main.nf" params(params)
include { multiqc_panel } from "./../../modules/multiqc/main.nf" params(params)

workflow PANEL_QC {

	take:
		bam
		panels
		targets

	main:
		panel_coverage(bam.combine(panels),targets.collect())
		multiqc_panel(panel_coverage.out[0])
	

	emit:
		qc = multiqc_panel.out

}
	
