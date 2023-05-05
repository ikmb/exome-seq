include { GATK_MUTECT2_PAIR } from "./../modules/gatk/mutect2_pair"
include { GATK_FILTER_MUTECT_CALLS } from "./../modules/gatk/filter_mutect_calls"
include { BCFTOOLS_VIEW } from "./../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_ANNOTATE } from "./../modules/bcftools/annotate"
include { GATK_LEARN_READ_ORIENTATION_MODEL } from "./../modules/gatk/learn_read_orientation_model"
include { GATK_GET_PILEUP_SUMMARIES } from "./../modules/gatk/get_pileup_summaries"
include { GATK_CALCULATE_CONTAMINATION } from "./../modules/gatk/calculate_contamination"

ch_versions 	= Channel.from([])
ch_vcfs 	= Channel.from([])

workflow GATK_MUTECT2_PAIRED {

	take:
		bams
		targets
		fasta
		dbsnp
		mutect_normals
		mutect_normals_tbi

	main:

		ch_bam_normal = bams.map { m,nb,ni,tb,ti -> 
			[ m,nb,ni ]
		}
		ch_bam_tumor = bams.map {  m,nb,ni,tb,ti ->
			[ m,tb,ti ]
		}
	
		GATK_MUTECT2_PAIR(
			bams.map { m,bn,bni,bt,bti ->
                            [
                                m,[bn,bt],[bni,bti]
                            ]
                        },
			targets.collect(),
			fasta.collect(),
			mutect_normals.collect(),
			mutect_normals_tbi.collect()
		)

		ch_versions = ch_versions.mix(GATK_MUTECT2_PAIR.out.versions)

		GATK_LEARN_READ_ORIENTATION_MODEL(
			GATK_MUTECT2_PAIR.out.f1r2
		)

		ch_versions = ch_versions.mix(GATK_LEARN_READ_ORIENTATION_MODEL.out.versions)

		GATK_GET_PILEUP_SUMMARIES(
			ch_bam_tumor,
			targets.collect(),
			fasta.collect()
		)

		ch_versions = ch_versions.mix(GATK_GET_PILEUP_SUMMARIES.out.versions)

		GATK_CALCULATE_CONTAMINATION(
			GATK_GET_PILEUP_SUMMARIES.out.table
		)

		ch_versions = ch_versions.mix(GATK_CALCULATE_CONTAMINATION.out.versions)

		ch_mutect = GATK_MUTECT2_PAIR.out.vcf.join(
			GATK_LEARN_READ_ORIENTATION_MODEL.out.model
		).join(GATK_CALCULATE_CONTAMINATION.out.table)

		GATK_FILTER_MUTECT_CALLS(
			ch_mutect,
			fasta.collect()
		)

		ch_versions = ch_versions.mix(GATK_FILTER_MUTECT_CALLS.out.versions)

		ch_vcfs = ch_vcfs.mix(GATK_FILTER_MUTECT_CALLS.out.vcf)

		BCFTOOLS_VIEW(ch_vcfs)

		ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions)

		BCFTOOLS_ANNOTATE_DBSNP(
			BCFTOOLS_VIEW.out.vcf.map { meta,v,t ->
                def s_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "MUTECT2" ]
                tuple(s_meta,v,t)
            },
			dbsnp.collect()
		)
		
		ch_versions = ch_versions.mix(BCFTOOLS_ANNOTATE_DBSNP.out.versions)

		BCFTOOLS_ANNOTATE(
			BCFTOOLS_ANNOTATE_DBSNP.out.vcf
        )

	emit:
	versions 	= ch_versions
	vcf 		= BCFTOOLS_ANNOTATE.out.vcf
		
}

