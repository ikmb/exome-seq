include { deepvariant ; merge_gvcfs } from "./../../modules/deepvariant/main.nf" params(params)
include { vcf_index ; vcf_add_dbsnp ; vcf_add_header ; vcf_filter_pass } from "./../../modules/vcf/main.nf" params(params)
include { merge_vcf }  from "./../../modules/vcf/main.nf"

workflow DV_VARIANT_CALLING {

	take:
		bam
		bed
		deepvariant_index

	main:
		deepvariant(
			bam,
			bed.collect(),
			deepvariant_index.collect()
		)
		vcf_index(deepvariant.out.vcf)
		vcf_filter_pass(vcf_index.out.vcf)
		vcf_add_dbsnp(vcf_filter_pass.out.vcf)
		vcf_add_header(vcf_add_dbsnp.out.vcf)
		// Joint calling with GLNexus - or simple merging
                if (params.joint_calling) {
                        merge_gvcfs(deepvariant.out.gvcf.collect(),bed)
                        merged_vcf = merge_gvcfs.out
                } else {
                        merge_vcf(deepvariant.out.vcf
				.map { m,v ->
					tuple(meta.caller,v)
				}
				.groupTuple().collect()
			)
                        merged_vcf = merge_vcf.out.vcf
                }

	emit:
		vcf = vcf_add_header.out
		sample_names = deepvariant.out.sample_name
		gvcf = deepvariant.out.gvcf
		vcf_multi = merged_vcf
}
