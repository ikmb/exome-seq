process GATK_SPLITINTERVALS {

	tag "${intervals}"

	label 'gatk'

	input:
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	path("results/*"), emit: intervals

	script:

	"""
		gatk SplitIntervals -R ${fasta} -L $intervals \
		--scatter-count 10 \
		-O results
	"""

}
