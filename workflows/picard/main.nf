include { hybrid_capture_metrics; multi_metrics; oxo_metrics } from "./../../modules/picard/main.nf" params(params)

workflow PICARD_METRICS {

	take:
		bam
		targets
		baits
	
	main:
		hybrid_capture_metrics(bam,targets.collect(),baits.collect())
		multi_metrics(bam,baits.collect())
		oxo_metrics(bam,targets.collect())

	emit:
		qc_reports = hybrid_capture_metrics.out[0].mix(multi_metrics.out[0],oxo_metrics.out[0])
}
