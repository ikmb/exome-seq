include { trim } from "./../../modules/trim/main.nf" params(params)
include { align } from "./../../modules/align/main.nf" params(params)
include { merge_multi_lane; bam_index ; dedup } from "./../../modules/samtools/main.nf" params(params)

workflow TRIM_AND_ALIGN {

	take:
		reads

	main:
		trim(reads)
		align( trim.out[0] )
		merge_multi_lane( align.out.groupTuple(by: [0,1]).filter{ i,s,b -> b.size() > 1 && b.size() < 1000 } )
		bam_index(merge_multi_lane.out.mix( align.out.groupTuple(by: [0,1]).filter{ i,s,b -> b.size() < 2 || b.size() > 1000 })  )
		dedup(bam_index.out[0])
		
	emit:
		bam = dedup.out[0]
		qc = trim.out[1]
		dedup_report = dedup.out[3]
}

