include { GUNZIP } from "./../modules/gunzip"
include { CNVKIT_BATCH } from "./../modules/cnvkit/batch"
include { CNVKIT_EXPORT } from "./../modules/cnvkit/export"

workflow CNVKIT {

	take:
		bam
		cnn_gz
	main:
		
		GUNZIP(cnn_gz)

		CNVKIT_BATCH(
			bam,
			GUNZIP.out.decompressed.collect()
		)

		CNVKIT_EXPORT(
			CNVKIT_BATCH.out.cns
		)
		
}
