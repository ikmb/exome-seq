include { GATK_MUTECT2_PAIR } from "./../modules/gatk/mutect2_pair"
include { GATK_FILTER_MUTECT_CALLS } from "./../modules/gatk/filter_mutect_calls"
include { BCFTOOLS_VIEW } from "./../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_ANNOTATE } from "./../modules/bcftools/annotate"
include { GATK_LEARN_READ_ORIENTATION_MODEL } from "./../modules/gatk/learn_read_orientation_model"
include { GATK_GET_PILEUP_SUMMARIES } from "./../modules/gatk/get_pileup_summaries"
include { GATK_CALCULATE_CONTAMINATION } from "./../modules/gatk/calculate_contamination"

workflow GATK_MUTECT2_PAIRED {

	take:
		bams
		targets
		fasta
		dbsnp
	
	main:

		ch_vcfs = Channel.from([])

		ch_bam_normal = bams.map { m,nb,ni,tb,ti -> 
			[ m,nb,ni ]
		}
		ch_bam_tumor = bams.map {  m,nb,ni,tb,ti ->
			[ m,tb,ti ]
		}
	
		GATK_MUTECT2_PAIR(
			bams,
			targets.collect(),
			fasta.collect()
		)

		GATK_LEARN_READ_ORIENTATION_MODEL(
			GATK_MUTECT2_PAIR.out.f1r2
		)

		GATK_GET_PILEUP_SUMMARIES(
			ch_bam_tumor,
			targets.collect()
		)

		GATK_CALCULATE_CONTAMINATION(
			GATK_GET_PILEUP_SUMMARIES.out.table
		)

		ch_mutect = GATK_MUTECT2_PAIR.out.vcf.join(
			GATK_LEARN_READ_ORIENTATION_MODEL.out.model
		).join(GATK_CALCULATE_CONTAMINATION.out.table)

		GATK_FILTER_MUTECT_CALLS(
			ch_mutect,
			fasta.collect()
		)

		ch_vcfs = ch_vcfs.mix(GATK_FILTER_MUTECT_CALLS.out.vcf)

		BCFTOOLS_VIEW(ch_vcfs)

		BCFTOOLS_ANNOTATE_DBSNP(
            BCFTOOLS_VIEW.out.vcf,
            dbsnp.collect()
        )
		
		BCFTOOLS_ANNOTATE(
			BCFTOOLS_ANNOTATE_DBSNP.out.vcf.map { meta,v,t ->
				def s_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "MUTECT2" ]
				tuple(s_meta,v,t)
			}
        )

	emit:

	vcf = BCFTOOLS_ANNOTATE.out.vcf
		
}

