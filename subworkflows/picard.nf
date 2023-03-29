include { HYBRID_CAPTURE_METRICS } from "./../modules/picard/collecthsmetrics"
include { MULTI_METRICS } from "./../modules/picard/collectmultiplemetrics"
include { OXO_METRICS } from "./../modules/picard/collectoxometrics"

ch_versions = Channel.from([])

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

		ch_versions = ch_versions.mix(HYBRID_CAPTURE_METRICS.out.versions)

		MULTI_METRICS(
			bam,
			baits.collect(),
			fasta.collect(),
			dbsnp.collect()

		)

		ch_versions = ch_versions.mix(MULTI_METRICS.out.versions)

		OXO_METRICS(
			bam,
			targets.collect(),
			fasta.collect(),
			dbsnp.collect()
		)

		ch_versions = ch_versions.mix(OXO_METRICS.out.versions)

	emit:
		versions = ch_versions
		qc_reports = HYBRID_CAPTURE_METRICS.out[0].mix(MULTI_METRICS.out[0],OXO_METRICS.out[0])
}
