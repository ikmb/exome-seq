include { GATK_BASERECALIBRATOR } from './../modules/gatk/baserecalibrator'
include { GATK_APPLYBQSR } from  './../modules/gatk/applybqsr'
include { GATK_GATHERBQSRREPORTS } from './../modules/gatk/gatherbqsrreports'
include { PICARD_SET_BAM_TAGS } from "./../modules/picard/setnmmdanduqtags"
include { GATK_SPLITINTERVALS } from './../modules/gatk/splitintervals'
include { SAMTOOLS_INDEX } from './../modules/samtools/index'

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
                fasta.collect()
        )

        // Parallelize BQSR recal computation
        GATK_SPLITINTERVALS.out.intervals.flatMap { i ->
                i.collect { file(it) }
        }.set { ch_intervals }

        // Add all relevant tags to BAM file
        PICARD_SET_BAM_TAGS(
                bam,
                fasta.collect()
        )
        SAMTOOLS_INDEX(
                PICARD_SET_BAM_TAGS.out.bam
        )

        // Recalibrate BAM base quality scores
        GATK_BASERECALIBRATOR(
                bam.combine(ch_intervals),
                fasta.collect(),
                known_snps.collect(),
                known_snps_tbi.collect(),
                known_indels.collect(),
                known_indels_tbi.collect()
        )

        recal_by_sample = GATK_BASERECALIBRATOR.out.report.groupTuple()

        GATK_GATHERBQSRREPORTS(
                recal_by_sample
        )

	// Apply recalibration
        GATK_APPLYBQSR(
                SAMTOOLS_INDEX.out.bam.join(GATK_GATHERBQSRREPORTS.out.report),
                intervals.collect(),
                fasta.collect()
        )
	
	emit:
	bam = GATK_APPLYBQSR.out.bam	

}
