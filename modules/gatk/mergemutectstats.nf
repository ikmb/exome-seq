process GATK_MERGEMUTECTSTATS {

    tag "${meta.patient_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2", mode: 'copy'    

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    input:
    tuple val(meta),path(all_stats)

    output:
    tuple val(meta),path(merged_stats), emit: stats
    path("versions.yml"), emit: versions

    script:
    merged_stats = meta.sample_id + "_mutect-merged.stats"

    """
    gatk --java-options "-Xmx4g" MergeMutectStats \
        --stats ${all_stats.join(' --stats ')} \
        -O $merged_stats
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
