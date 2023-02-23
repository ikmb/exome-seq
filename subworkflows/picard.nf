include { HYBRID_CAPTURE_METRICS } from "./../modules/picard/collecthsmetrics"
include { MULTI_METRICS } from "./../modules/picard/collectmultiplemetrics"
include { OXO_METRICS } from "./../modules/picard/collectoxometrics"

workflow PICARD_METRICS {

	take:
		bam
		targets
		baits
		fasta
		dbsnp
	
	main:
		HYBRID_CAPTURE_METRICS(
			bam,targets.collect(),
			baits.collect(),
			fasta.collect()
		)
		MULTI_METRICS(
			bam,
			baits.collect(),
			fasta.collect(),
			dbsnp.collect()

		)
		OXO_METRICS(
			bam,
			targets.collect(),
			fasta.collect(),
			dbsnp.collect()
		)

	emit:
		qc_reports = HYBRID_CAPTURE_METRICS.out[0].mix(MULTI_METRICS.out[0],OXO_METRICS.out[0])
}
