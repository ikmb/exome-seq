include { FASTP as TRIM } from "./../modules/fastp"
include { BWA_MEM } from "./../modules/bwa/mem"
include { BWA2_MEM } from "./../modules/bwa2/mem"
include { DRAGMAP_ALIGN } from "./../modules/dragmap/align"
include { SAMTOOLS_MERGE as MERGE_MULTI_LANE } from "./../modules/samtools/merge" 
include { SAMTOOLS_INDEX as BAM_INDEX; SAMTOOLS_INDEX as BAM_INDEX_FILTERED } from "./../modules/samtools/index"
include { SAMTOOLS_MARKDUP as DEDUP } from "./../modules/samtools/markdup"
include { SAMTOOLS_AMPLICONCLIP as AMPLICON_CLIP } from "./../modules/samtools/ampliconclip"

ch_align_log = Channel.from([])

workflow TRIM_AND_ALIGN {

	take:
		samplesheet
		amplicon_bed
		genome_index
	main:

		Channel.fromPath(samplesheet)
		.splitCsv ( header: true, sep: ';')
		.map { create_fastq_channel(it) }
		.set {reads }

		TRIM(
			reads
        )

		ch_aligned_bams = Channel.from([])

		// run Dragen aligner or BWA/BWA2
		if (params.aligner == "dragmap") {
			DRAGMAP_ALIGN(
				TRIM.out.reads,
				genome_index
			)
			ch_aligned_bams = ch_aligned_bams.mix(DRAGMAP_ALIGN.out.bam)
			ch_sample_names = DRAGMAP_ALIGN.out.sample_name
			ch_align_log = ch_align_log.mix(DRAGMAP_ALIGN.out.log)
		} else if (params.aligner == "bwa") {
			
			BWA_MEM( 
				TRIM.out.reads, 
				genome_index
			)
			ch_aligned_bams = ch_aligned_bams.mix(BWA_MEM.out.bam)
			ch_sample_names = BWA_MEM.out.sample_name
		} else if (params.aligner == "bwa2") {
			BWA2_MEM(
                                TRIM.out.reads,
                                genome_index
                        )
                        ch_aligned_bams = ch_aligned_bams.mix(BWA2_MEM.out.bam)
                        ch_sample_names = BWA2_MEM.out.sample_name
		}

		bam_mapped = ch_aligned_bams.map { meta, bam ->
                        new_meta = [:]
			new_meta.patient_id = meta.patient_id
			new_meta.sample_id = meta.sample_id
			def groupKey = meta.sample_id
			tuple( groupKey, new_meta, bam)
		}.groupTuple(by: [0,1]).map { g ,new_meta ,bam -> [ new_meta, bam ] }
			
		bam_mapped.branch {
			single:   it[1].size() == 1
			multiple: it[1].size() > 1
		}.set { bam_to_merge }

		MERGE_MULTI_LANE( bam_to_merge.multiple )
		BAM_INDEX(MERGE_MULTI_LANE.out.bam.mix( bam_to_merge.single ))

		ch_report = Channel.from([])
		if (params.amplicon_bed) {
			//AMPLICON_CLIP(
			//	BAM_INDEX.out.bam,
			//	amplicon_bed.collect()
			//)
			//ch_final_bam = AMPLICON_CLIP.out.bam
			ch_final_bam = BAM_INDEX.out.bam
		} else {
			DEDUP(BAM_INDEX.out.bam)
			ch_final_bam = DEDUP.out.bam
			ch_report = ch_report.mix(DEDUP.out.report)
		}
		
	emit:
		bam_nodedup = BAM_INDEX.out.bam
		bam = ch_final_bam
		qc = TRIM.out.json
		logs = ch_align_log
		dedup_report = ch_report
		sample_names =ch_sample_names.unique()
		metas = MERGE_MULTI_LANE.out.meta_data
}

def create_fastq_channel(LinkedHashMap row) {

    // IndivID;SampleID;libraryID;rgID;rgPU;platform;platform_model;Center;Date;R1;R2

    def meta = [:]
    meta.patient_id = row.IndivID
    meta.sample_id = row.SampleID
    meta.library_id = row.libraryID
    meta.readgroup_id = row.rgID
    meta.center = row.Center
    meta.date = row.Date
    meta.platform_unit = row.rgPU

    def array = []
    array = [ meta, file(row.R1), file(row.R2) ]

    return array
}

