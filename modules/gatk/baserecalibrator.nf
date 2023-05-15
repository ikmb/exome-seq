process GATK_BASERECALIBRATOR {

	tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'medium_serial'

	input:
	tuple val(meta),path(bam),path(bai),path(intervals)
	tuple path(fasta),path(fai),path(dict)
	path(snps)
	path(snps_tbi)
	path(indels)
	path(indels_tbi)

	output:
	tuple val(meta),path(report), emit: report
	path("versions.yml"), emit: versions

	script:
	report = bam.getBaseName() + "-" + intervals.getBaseName() + "-recal.txt"	

    """
    gatk  --java-options "-Xmx${task.memory.giga}g" BaseRecalibrator -R ${fasta}  -I $bam -O $report \
        -L $intervals --use-original-qualities --known-sites ${snps.join(' --known-sites ')} \
        --known-sites ${indels.join(' --known-sites ')} \
        -ip $params.interval_padding

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

    """
}
