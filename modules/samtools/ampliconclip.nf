process SAMTOOLS_AMPLICONCLIP {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'	

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/", mode: 'copy'

    input:
    tuple val(meta),path(bam),path(bai)
    path(bed)

    output:
    tuple val(meta),path(bam_masked),path(bam_masked_bai), emit: bam
    path("versions.yml"), emit: versions

    script:
    bam_masked = bam.getBaseName() + "-amplicon_clipped.bam"
    bam_masked_bai = bam_masked + ".bai"

    """
    samtools ampliconclip -b $bed $bam | samtools sort -o $bam_masked
    samtools index $bam_masked

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}


