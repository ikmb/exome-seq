include { interval_to_bed } from "./../../modules/picard/main.nf" params(params)
include { bed_to_bedgz } from "./../../modules/bed/main.nf" params(params)

workflow CONVERT_BED {
	
	take:
		intervals

	main:
		interval_to_bed(intervals)
		bed_to_bedgz(interval_to_bed.out)

	emit:
		bed = interval_to_bed.out
		bed_gz = bed_to_bedgz.out

}
