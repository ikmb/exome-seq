process GATK_MAKESITESONLYVCF {

	label 'gatk'

	input:
	tuple val(meta),path(vcf_file),path(tbi)

	output:
	tuple val(meta),path(vcf_sites_only),path(vcf_sites_only_tbi), emit: vcf

	script:
	vcf_sites_only = vcf_file.getSimpleName() + "-sites_only.vcf.gz"
	vcf_sites_only_tbi = vcf_sites_only + ".tbi"

	"""
		gatk MakeSitesOnlyVcf -I $vcf_file -O $vcf_sites_only
		gatk IndexFeatureFile -I $vcf_sites_only
	"""
}
