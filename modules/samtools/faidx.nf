process SAMTOOLS_FAIDX {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    publishDir "${params.outdir}/${assembly}", mode: 'copy'

    tag "${fasta}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta),path(fai), emit: fai
    path("versions.yml"), emit: versions

    script:
    assembly = meta.assembly
    fai = fasta + ".fai"

    """
    samtools faidx $fasta > $fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

