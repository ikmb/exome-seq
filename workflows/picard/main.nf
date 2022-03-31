include { HYBRID_CAPTURE_METRICS ; MULTI_METRICS ; OXO_METRICS } from "./../../modules/picard/main.nf" params(params)

workflow PICARD_METRICS {

	take:
		bam
		targets
		baits
	
	main:
		HYBRID_CAPTURE_METRICS(bam,targets.collect(),baits.collect())
		MULTI_METRICS(bam,baits.collect())
		OXO_METRICS(bam,targets.collect())

	emit:
		qc_reports = HYBRID_CAPTURE_METRICS.out[0].mix(MULTI_METRICS.out[0],OXO_METRICS.out[0])
}
