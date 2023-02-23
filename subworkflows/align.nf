include { FASTP as TRIM } from "./../modules/fastp.nf"
include { BWA_MEM as ALIGN } from "./../modules/bwa/mem.nf"
include { SAMTOOLS_MERGE as MERGE_MULTI_LANE } from "./../modules/samtools/merge.nf" 
include { SAMTOOLS_INDEX as BAM_INDEX; SAMTOOLS_INDEX as BAM_INDEX_FILTERED } from "./../modules/samtools/index.nf"
include { SAMTOOLS_MARKDUP as DEDUP } from "./../modules/samtools/markdup.nf"
include { SAMTOOLS_AMPLICONCLIP as AMPLICON_CLIP } from "./../modules/samtools/ampliconclip.nf"

workflow TRIM_AND_ALIGN {

	take:
		samplesheet
		amplicon_bed
		bwa_index
	main:

		Channel.fromPath(samplesheet)
		.splitCsv ( header: true, sep: ';')
		.map { create_fastq_channel(it) }
		.set {reads }

		TRIM(
			reads
                )
		ALIGN( 
			TRIM.out.reads, 
			bwa_index
		)
		bam_mapped = ALIGN.out.bam.map { meta, bam ->
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
		dedup_report = ch_report
		sample_names = ALIGN.out.sample_name.unique()
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

