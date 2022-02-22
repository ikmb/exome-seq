include { trim } from "./../../modules/trim/main.nf" params(params)
include { align } from "./../../modules/align/main.nf" params(params)
include { merge_multi_lane; bam_index ; dedup } from "./../../modules/samtools/main.nf" params(params)

workflow TRIM_AND_ALIGN {

	take:
		samplesheet

	main:

		Channel.fromPath(samplesheet)
		.splitCsv ( header: true, sep: ';')
		.map { create_fastq_channel(it) }
		.set {reads }

		trim(
			reads
                )
		align( trim.out.reads )
		bam_mapped = align.out.bam.map { meta, bam ->
			tuple( meta.sample_id, meta, bam)
		}.groupTuple().map { mid,meta,bam -> [ meta, bam ] }
			
		bam_mapped.branch {
		        single:   it[1].size() == 1
		        multiple: it[1].size() > 1
	        }.set { bam_to_merge }

		merge_multi_lane( bam_to_merge.multiple )
		bam_index(merge_multi_lane.out.bam.mix( bam_to_merge.single ))
		dedup(bam_index.out.bam)
		
	emit:
		bam = dedup.out.bam
		qc = trim.out.json
		dedup_report = dedup.out.report
		sample_names = align.out.sample_name.unique()
}

def create_fastq_channel(LinkedHashMap row) {

    // IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2

    def meta = [:]
    meta.patient_id = row.indivID
    meta.sample_id = row.sampleID
    meta.library_id = row.libraryID
    meta.readgroup_id = row.rgID
    meta.center = row.Center
    meta.date = row.Date
    meta.platform_unit = row.rgPU

    def array = []
    array = [ meta, file(row.R1), file(row.R2) ]

    return array
}

