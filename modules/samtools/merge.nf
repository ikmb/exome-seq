process SAMTOOLS_MERGE {

    container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta), path(aligned_bam_list)

    output:
    tuple val(meta),path(merged_bam), emit: bam
    val(meta), emit: meta_data
    path("versions.yml"), emit: versions

    script:
    merged_bam = meta.patient_id + "_" + meta.sample_id + "-merged.bam"
    merged_bam_index = merged_bam + ".bai"
    
    """
    samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

