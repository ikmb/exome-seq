process BCFTOOLS_STATS {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Stats", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)

    output:
    path(vcf_stats), emit: stats
    path("versions.yml"), emit: versions

    script:
    vcf_stats = vcf.getBaseName() + ".stats"

    """
    bcftools stats $vcf > $vcf_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}

