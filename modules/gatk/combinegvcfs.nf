process GATK_COMBINEGVCFS {

	label 'gatk'

	input:
	path(gvcfs)
	path(tbis)
	path(intervals)
	path(fasta)

	output:
	tuple path(merged_gvcf),path(merged_gvcf_tbi), emit: gvcf

	script:
	merged_gvcf = "GATK." + params.run_name + ".merged.g.vcf.gz"
	merged_gvcf_tbi = merged_gvcf + ".tbi"

	"""
		gatk CombineGVCFs -R $fasta --variant ${gvcfs.join(' --variant ')} \
			-O $merged_gvcf -OVI true -L $intervals -ip $params.interval_padding
	"""

}
