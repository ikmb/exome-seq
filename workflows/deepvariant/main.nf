include { DEEPVARIANT ; MERGE_GVCFS } from "./../../modules/deepvariant/main.nf" params(params)
include { VCF_INDEX ; VCF_ADD_DBSNP ; VCF_ADD_HEADER ;  VCF_FILTER_PASS } from "./../../modules/vcf/main.nf" params(params)
include { MERGE_VCF }  from "./../../modules/vcf/main.nf"
include { WHATSHAP; WHATSHAP_SINGLE } from "./../../modules/whatshap/main.nf"

workflow DV_VARIANT_CALLING {

	take:
		bam
		bed
		deepvariant_index

	main:
                merged_vcf = Channel.empty()
		ch_phased_multi = Channel.empty()

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

		// Phase single VCF
		WHATSHAP_SINGLE(
			VCF_ADD_HEADER.out.vcf.map {m,v,t ->
				new_key = m.sample_id
				tuple(new_key,m,v,t)
			}.join(
				bam.map { m,b,i ->
					new_b_key = m.sample_id
					tuple(new_b_key,b,i)
				}
			).map { n,m,v,t,b,i -> [ m,v,t,b,i ] }
		)

		// Fiddly work-around to determine whether we have 1 or multiple vcfs. No merging when n=1
		VCF_FILTER_PASS.out.vcf.map { m,v,t ->
			new_meta = [ id: "all", sample_id: "UNDEFINED", patient_id: "UNDEFINED", variantcaller: "DEEPVARIANT" ]
			tuple(new_meta,v,t)
		}.set { ch_variants_pass }

		ch_variants_pass
		.groupTuple()
		.branch { m,v,t ->
			single: v.size() == 1
			multi: v.size() > 1
		}.set { ch_grouped_vcfs }
		
		// Joint calling with GLNexus - or simple merging
                if (params.joint_calling) {
	                MERGE_GVCFS(DEEPVARIANT.out.gvcf.collect(),bed)
        	        joint_vcf = MERGE_GVCFS.out
                        def j_meta = [ id: "JointCalling", sample_id: "GLNEXUS_DEEPVARIANT", patient_id: "MergedCallset", variantcaller: "DEEPVARIANT" ]
                        merged_vcf = merged_vcf.mix(
				joint_vcf.map { v,t -> tuple(j_meta,v,t) }
			)	
			WHATSHAP(
				merged_vcf,
				bam.map{m,b,i -> tuple(b,i) }.collect()
			)
			ch_phased_multi = ch_phased_multi.mix(WHATSHAP.out.vcf)
	        } else {
       	                MERGE_VCF(
  				ch_grouped_vcfs.multi.map { m,v,t -> [ [ id: "Deepvariant", sample_id: "Bcftools", patient_id: "MergedCallset", variantcaller: "DEEPVARIANT"],v,t] }
			)
               	        merged_vcf = merged_vcf.mix(MERGE_VCF.out.vcf)
	        }

	emit:
		vcf = VCF_ADD_HEADER.out.vcf
		vcf_phased_single = WHATSHAP_SINGLE.out.vcf
		sample_names = DEEPVARIANT.out.sample_name
		gvcf = DEEPVARIANT.out.gvcf
		vcf_multi = merged_vcf
		vcf_phased_multi = ch_phased_multi
}
