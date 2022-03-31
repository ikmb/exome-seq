include { INTERVAL_TO_BED } from "./../../modules/picard/main.nf" params(params)
include { BED_TO_BEDGZ } from "./../../modules/bed/main.nf" params(params)

workflow CONVERT_BED {
	
	take:
		intervals

	main:
		INTERVAL_TO_BED(intervals)
		BED_TO_BEDGZ(INTERVAL_TO_BED.out)

	emit:
		bed = INTERVAL_TO_BED.out
		bed_gz = BED_TO_BEDGZ.out

}
