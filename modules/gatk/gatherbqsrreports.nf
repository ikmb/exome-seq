process GATK_GATHERBQSRREPORTS {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'short_serial'

    input:
    tuple val(meta),path(reports)

    output:
    tuple val(meta),path(merged_report), emit: report
    path("versions.yml"), emit: versions

    script:
    merged_report = meta.patient_id + "_" + meta.sample_id + "-recal.txt"

    """
    gatk GatherBQSRReports --input ${reports.join(' --input ')} --output $merged_report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
