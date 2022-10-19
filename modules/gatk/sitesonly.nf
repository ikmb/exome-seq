process GATK_MAKESITESONLYVCF {

	label 'gatk'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_sites_only),path(vcf_sites_only_tbi), emit: vcf

	script:
	vcf_sites_only = vcf.getBaseName() + ".sites_only.vcf.gz"
	vcf_sites_only_tbi = vcf_sites_only + ".tbi"

	"""
		gatk MakeSitesOnlyVcf -I $vcf -O $vcf_sites_only
		gatk IndexFeatureFile -I $vcf
	"""
}
