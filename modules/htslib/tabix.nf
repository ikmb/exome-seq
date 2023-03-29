process TABIX {

    container 'quay.io/biocontainers/htslib:1.16--h6bc39ce_0'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(vcf),path(tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:

    tbi = vcf + ".tbi"

    """
    tabix $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}

