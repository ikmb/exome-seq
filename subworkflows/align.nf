include { FASTP as TRIM } from "./../modules/fastp"
include { BWA_MEM } from "./../modules/bwa/mem"
include { BWA2_MEM } from "./../modules/bwa2/mem"
include { DRAGMAP_ALIGN } from "./../modules/dragmap/align"
include { SAMTOOLS_MERGE as MERGE_MULTI_LANE } from "./../modules/samtools/merge" 
include { SAMTOOLS_INDEX as BAM_INDEX; SAMTOOLS_INDEX as BAM_INDEX_FILTERED } from "./../modules/samtools/index"
include { SAMTOOLS_MARKDUP as DEDUP } from "./../modules/samtools/markdup"
include { SAMTOOLS_AMPLICONCLIP } from "./../modules/samtools/ampliconclip"
include { SAMTOOLS_AMPLICONSTATS } from "./../modules/samtools/ampliconstats"
include { SAMTOOLS_BAM2CRAM } from "./../modules/samtools/bam2cram"
include { BAM_MD5 } from "./../modules/bam_md5"
include { VALIDATE_AMPLICONBED } from "./../modules/validate_ampliconbed.nf"

ch_align_log = Channel.from([])
ch_versions = Channel.from([])

workflow TRIM_AND_ALIGN {

    take:
        samplesheet
        amplicon_bed
        genome_index
        fasta

    main:

        samplesheet
        .splitCsv ( header: true, sep: ';')
        .map { create_fastq_channel(it) }
        .set {reads }

        TRIM(
            reads
        )

        ch_versions = ch_versions.mix(TRIM.out.versions)

        ch_aligned_bams = Channel.from([])

        // run Dragen aligner or BWA/BWA2
        if (params.aligner == "dragmap") {
            DRAGMAP_ALIGN(
                TRIM.out.reads,
                genome_index
            )
            ch_aligned_bams = ch_aligned_bams.mix(DRAGMAP_ALIGN.out.bam)
            ch_sample_names = DRAGMAP_ALIGN.out.sample_name
            ch_align_log    = ch_align_log.mix(DRAGMAP_ALIGN.out.log)

            ch_versions = ch_versions.mix(DRAGMAP_ALIGN.out.versions)

        } else if (params.aligner == "bwa") {
            
            BWA_MEM( 
                TRIM.out.reads, 
                genome_index
            )
            ch_aligned_bams = ch_aligned_bams.mix(BWA_MEM.out.bam)
            ch_sample_names = BWA_MEM.out.sample_name

            ch_versions = ch_versions.mix(BWA_MEM.out.versions)

        } else if (params.aligner == "bwa2") {
            BWA2_MEM(
                TRIM.out.reads,
                genome_index
            )
            ch_aligned_bams = ch_aligned_bams.mix(BWA2_MEM.out.bam)
            ch_sample_names = BWA2_MEM.out.sample_name

            ch_versions = ch_versions.mix(BWA2_MEM.out.versions)
        }

        bam_mapped = ch_aligned_bams.map { meta, bam ->
            new_meta = [:]
            new_meta.patient_id = meta.patient_id
            new_meta.sample_id = meta.sample_id
            new_meta.status = meta.status
            def groupKey = meta.sample_id
            tuple( groupKey, new_meta, bam)
        }.groupTuple(by: [0,1]).map { g ,new_meta ,bam -> [ new_meta, bam ] }
            
        bam_mapped.branch {
            single:   it[1].size() == 1
            multiple: it[1].size() > 1
        }.set { bam_to_merge }

        MERGE_MULTI_LANE( bam_to_merge.multiple )

        ch_versions = ch_versions.mix(MERGE_MULTI_LANE.out.versions)

        BAM_INDEX(MERGE_MULTI_LANE.out.bam.mix( bam_to_merge.single ))

        ch_report = Channel.from([])

        // Data is from amplicon sequencing, remove amplicon artifacts
        if (params.amplicon_bed) {

            VALIDATE_AMPLICONBED(
                amplicon_bed
            )

            SAMTOOLS_AMPLICONCLIP(
                BAM_INDEX.out.bam,
                VALIDATE_AMPLICONBED.out.bed.collect()
            )

            ch_final_bam = SAMTOOLS_AMPLICONCLIP.out.bam

            SAMTOOLS_AMPLICONSTATS(
                ch_final_bam,
                VALIDATE_AMPLICONBED.out.bed.collect()
            )
            
            //ch_final_bam = BAM_INDEX.out.bam

        } else {

            DEDUP(
                BAM_INDEX.out.bam,
                fasta
            )
            
            ch_final_bam     = DEDUP.out.bam
            ch_report        = ch_report.mix(DEDUP.out.report)
            ch_versions      = ch_versions.mix(DEDUP.out.versions)
        }

        //convert bam to cram and stage out
        SAMTOOLS_BAM2CRAM(
            ch_final_bam,
            fasta.collect()
        )
        
        // This creates the md5 sum AND stages the BAM file into the result folder. 
        // Else we need some more detailed logic to determine when a BAM file can be staged (amplicon vs other)
        BAM_MD5(
            SAMTOOLS_BAM2CRAM.out.cram
        )

    emit:
        bam_nodedup     = BAM_INDEX.out.bam
        bam             = ch_final_bam
        cram            = SAMTOOLS_BAM2CRAM.out.cram
        qc              = TRIM.out.json
        logs            = ch_align_log
        dedup_report    = ch_report
        sample_names    = ch_sample_names.unique()
        metas           = MERGE_MULTI_LANE.out.meta_data
        versions        = ch_versions
}

def create_fastq_channel(LinkedHashMap row) {

    // patient;sample;status;library;readgroup;platform_unit;center;date;R1;R2

    def meta = [:]
    
    meta.patient_id = row.patient
    meta.status = row.status.toInteger()
    meta.sample_id = row.sample
    meta.library_id = row.library
    meta.readgroup_id = row.readgroup
    meta.center = row.center
    meta.date = row.date
    meta.platform_unit = row.platform_unit

    def array = []
    array = [ meta, file(row.R1), file(row.R2) ]

    return array
}

