process SAMTOOLS_AMPLICONSTATS {

    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'	

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Processing/Ampliconstats", mode: 'copy'

    input:
    tuple val(meta),path(bam),path(bai)
    path(bed)

    output:
    tuple val(meta),path(amplicon_stats), emit: stats
    path("versions.yml"), emit: versions

    script:
    amplicon_stats = bam.getBaseName() + "-amplicon_stats.txt"

    """
    samtools ampliconstats $bed $bam > $amplicon_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}


