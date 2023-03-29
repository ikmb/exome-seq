process GATK_SPLITINTERVALS {

	tag "${intervals}"

	label 'gatk'

	input:
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	path("results/*"), emit: intervals
	path("versions.yml"), emit: versions

	script:

    """
    gatk SplitIntervals -R ${fasta} -L $intervals \
        --scatter-count 10 \
        -O results
		
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
	"""

}
