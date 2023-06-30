process SAMTOOLS_AMPLICONCLIP {

    container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'	

    label 'medium_serial' 

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/", mode: 'copy'

    input:
    tuple val(meta),path(bam_raw),path(bai_raw)
    path(bed)

    output:
    tuple val(meta),path(bam_masked),path(bam_masked_bai), emit: bam
    path("versions.yml"), emit: versions

    script:
    bam_masked = bam_raw.getBaseName() + "-amplicon_clipped.bam"
    bam_masked_bai = bam_masked + ".bai"

    """
    samtools ampliconclip -b $bed $bam_raw | samtools sort -o $bam_masked
    samtools index $bam_masked

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

}


