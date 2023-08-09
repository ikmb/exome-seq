process SAMTOOLS_BAM2CRAM {

    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta), path(cram),path(crai), emit: cram
    path("versions.yml"), emit: versions

    script:
    cram = bam.getBaseName() + ".cram"
    crai = cram + ".crai"

    """
    samtools view -O CRAM -o $cram -T $fasta $bam
    samtools index $cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}

