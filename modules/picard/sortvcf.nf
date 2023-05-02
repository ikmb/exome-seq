process PICARD_SORTVCF {

    label 'picard'

    label 'short_serial'

    input:
    tuple val(meta),path(vcf),path(tbi)

    output:
    tuple val(meta),path(vcf_sorted),path(tbi_sorted), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_sorted = vcf.getSimpleName() + "-sorted.vcf.gz"
    tbi_sorted = vcf_sorted + ".tbi"

    """
    picard SortVcf I=$vcf O=$vcf_sorted CREATE_INDEX=true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS

    """

}

