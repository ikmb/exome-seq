include { deepvariant } from "./../../modules/deepvariant/main.nf" params(params)
include { vcf_index ; vcf_add_dbsnp ; vcf_add_header ; vcf_filter_pass } from "./../../modules/vcf/main.nf" params(params)

workflow VARIANT_CALLING {

	take:
		bam
		bed
		deepvariant_index

	main:
		deepvariant(bam,bed.collect(),deepvariant_index.collect())
		vcf_index(deepvariant.out[1])
		vcf_filter_pass(vcf_index.out)
		vcf_add_dbsnp(vcf_filter_pass.out)
		vcf_add_header(vcf_add_dbsnp.out)
		
	emit:
		vcf = vcf_add_header.out
		sample_names = deepvariant.out[3]
		gvcf = deepvariant.out[0]
		vcf_only = deepvariant.out[2]
}
