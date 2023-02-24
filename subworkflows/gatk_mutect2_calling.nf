include { GATK_MUTECT2 } from "./../modules/gatk/mutect2"
include { BCFTOOLS_VIEW } from "./../modules/bcftools/view"
include { BCFTOOLS_ANNOTATE_DBSNP } from "./../modules/bcftools/annotate_dbsnp"
include { BCFTOOLS_ANNOTATE } from "./../modules/bcftools/annotate"

workflow GATK_MUTECT2_CALLING {

	take:
		bam
		targets
		fasta
		dbsnp
	
	main:
		
		GATK_MUTECT2(
			bam,
			targets.collect(),
			fasta.collect()
		)

		BCFTOOLS_VIEW(GATK_MUTECT2.out.vcf)

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

