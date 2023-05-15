include { GATK_BASERECALIBRATOR } from './../modules/gatk/baserecalibrator'
include { GATK_APPLYBQSR } from  './../modules/gatk/applybqsr'
include { GATK_GATHERBQSRREPORTS } from './../modules/gatk/gatherbqsrreports'
include { PICARD_SET_BAM_TAGS } from "./../modules/picard/setnmmdanduqtags"
include { GATK_SPLITINTERVALS } from './../modules/gatk/splitintervals'
include { SAMTOOLS_INDEX } from './../modules/samtools/index'

ch_versions = Channel.from([])

workflow GATK_BAM_RECAL {

	take:
	bam
	intervals
	fasta
	known_snps
	known_snps_tbi
	known_indels
	known_indels_tbi
	
	main:

        GATK_SPLITINTERVALS(
                intervals,
                fasta
        )

	ch_versions = ch_versions.mix(GATK_SPLITINTERVALS.out.versions)

        // Parallelize BQSR recal computation
        GATK_SPLITINTERVALS.out.intervals.flatMap { i ->
                i.collect { file(it) }
        }.set { ch_intervals }

        // Add all relevant tags to BAM file
        PICARD_SET_BAM_TAGS(
                bam,
                fasta
        )

	ch_versions = ch_versions.mix(PICARD_SET_BAM_TAGS.out.versions)

        SAMTOOLS_INDEX(
                PICARD_SET_BAM_TAGS.out.bam
        )

	ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

        // Recalibrate BAM base quality scores
        GATK_BASERECALIBRATOR(
                bam.combine(ch_intervals),
                fasta,
                known_snps,
                known_snps_tbi,
                known_indels,
                known_indels_tbi
        )

	ch_versions = ch_versions.mix(GATK_BASERECALIBRATOR.out.versions)

        recal_by_sample = GATK_BASERECALIBRATOR.out.report.groupTuple()

        GATK_GATHERBQSRREPORTS(
                recal_by_sample
        )

	ch_versions = ch_versions.mix(GATK_GATHERBQSRREPORTS.out.versions)

	// Apply recalibration
        GATK_APPLYBQSR(
                SAMTOOLS_INDEX.out.bam.join(GATK_GATHERBQSRREPORTS.out.report),
                intervals,
                fasta
        )

	ch_versions = ch_versions.mix(GATK_APPLYBQSR.out.versions)
	
	emit:
	bam = GATK_APPLYBQSR.out.bam	
	versions = ch_versions
}
