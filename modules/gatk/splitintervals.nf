process GATK_SPLITINTERVALS {

	tag "${intervals}"

	label 'gatk'

	input:
	path(intervals)

	output:
	path("results/*"), emit: intervals

	script:

	"""
		gatk SplitIntervals -R $params.fasta -L $intervals \
		--scatter-count 10 \
		-O results
	"""

}
