include { DEEPVARIANT ; MERGE_GVCFS } from "./../../modules/deepvariant/main.nf" params(params)
include { VCF_INDEX ; VCF_ADD_DBSNP ; VCF_ADD_HEADER ;  VCF_FILTER_PASS } from "./../../modules/vcf/main.nf" params(params)
include { MERGE_VCF }  from "./../../modules/vcf/main.nf"

workflow DV_VARIANT_CALLING {

	take:
		bam
		bed
		deepvariant_index

	main:
		
		DEEPVARIANT(
			bam,
			bed.collect(),
			deepvariant_index.collect()
		)

		VCF_INDEX(DEEPVARIANT.out.vcf)
		VCF_FILTER_PASS(VCF_INDEX.out.vcf)
		VCF_ADD_DBSNP(VCF_FILTER_PASS.out.vcf)
		VCF_ADD_HEADER(VCF_ADD_DBSNP.out.vcf.map { meta,v,t ->
			def s_meta = [ id: meta.id, sample_id: meta.sample_id, patient_id: meta.patient_id, variantcaller: "DEEPVARIANT" ]
			tuple(s_meta,v,t)
			}
		)

		// Joint calling with GLNexus - or simple merging
                if (params.joint_calling) {
	                MERGE_GVCFS(DEEPVARIANT.out.gvcf.collect(),bed)
        	        merged_vcf = MERGE_GVCFS.out
                        def j_meta = [ id: "JointCalling", sample_id: "GLNEXUS_DEEPVARIANT", patient_id: "MergedCallset", variantcaller: "DEEPVARIANT" ]
                        merged_vcf = merged_vcf.map { v,t -> [ j_meta,v,t]}
	        } else {
       	                MERGE_VCF(
				VCF_FILTER_PASS.out.vcf.map { m,v,t -> [ [ id: "Deepvariant", sample_id: "Bcftools", patient_id: "MergedCallset", variantcaller: "DEEPVARIANT"],v,t ] }.groupTuple()
			)
               	        merged_vcf = MERGE_VCF.out.vcf
	        }

	emit:
		vcf = VCF_ADD_HEADER.out.vcf
		sample_names = DEEPVARIANT.out.sample_name
		gvcf = DEEPVARIANT.out.gvcf
		vcf_multi = merged_vcf
}
