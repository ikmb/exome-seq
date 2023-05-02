process SAMTOOLS_DICT {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    publishDir "${params.outdir}/${assembly}", mode: 'copy'

    tag "${fasta}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta),path(dict), emit: dict
    path("versions.yml"), emit: versions

    script:
    assembly = meta.assembly
    dict = fasta.getBaseName() + ".dict"

    """
    samtools dict $fasta > $dict

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

