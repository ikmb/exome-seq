process concat {

	label 'bcftools'

	input:
	path(vcfs)

	output:
	path(merged_vcf), emit: vcf

	script:

	merged_vcf = "merged_callset." + params.run_name + ".vcf.gz"

	"""
		bcftools concat -a -D -o $merged_vcf -O z $vcfs
	"""

}
