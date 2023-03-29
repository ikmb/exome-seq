include { INTERVAL_TO_BED } from "./../modules/picard/intervallisttools" 
include { BGZIP_INDEX as BED_TO_BEDGZ } from "./../modules/htslib/bgzip_index"

workflow CONVERT_BED {
	
	take:
		intervals

	main:
		INTERVAL_TO_BED(intervals)
		BED_TO_BEDGZ(INTERVAL_TO_BED.out.bed)

	emit:
		bed = INTERVAL_TO_BED.out.bed
		bed_gz = BED_TO_BEDGZ.out.bedgz

}
