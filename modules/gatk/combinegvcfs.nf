process GATK_COMBINEGVCFS {

	tag "ALL"
	
	container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'medium_serial'

	input:
	path(gvcfs)
	path(tbis)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple path(merged_gvcf),path(merged_gvcf_tbi), emit: gvcf
	path("versions.yml"), emit: versions

	script:
	merged_gvcf = "GATK_" + params.run_name + "-merged.g.vcf.gz"
	merged_gvcf_tbi = merged_gvcf + ".tbi"

    """
    gatk CombineGVCFs -R $fasta --variant ${gvcfs.join(' --variant ')} \
        -O $merged_gvcf -OVI true -L $intervals -ip $params.interval_padding

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
	"""

}
