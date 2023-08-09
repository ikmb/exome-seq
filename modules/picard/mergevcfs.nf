process PICARD_MERGEVCFS {

    label 'picard'

    label 'short_serial'

    input:
    tuple val(meta),path(vcfs),path(tbis)
    output:
    tuple val(meta),path(vcf_merged),path(tbi_merged), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_merged = meta.sample_id + "-merged.vcf.gz"
    tbi_merged = vcf_merged + ".tbi"

    """
    picard MergeVcfs I=${vcfs.join(' -I ')} O=$vcf_merged CREATE_INDEX=true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS

    """

}

