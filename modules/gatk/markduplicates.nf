process GATK_MARK_DUPLICATES {

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    tag "${meta.patient_id}|${meta.sample_id}"

    input:
    tuple val(meta),path(bam),path(bai)

    output:
    tuple val(meta),path(bam_md),path(bai_md), emit: bam
    tuple val(meta),path(stats), emit: stats
    path("versions.yml"), emit: versions

    script:

    bam_md = bam.getBaseName() + "-dedup.bam"
    bai_md = bam.getBaseName() + "-dedup.bam.bai"
    metrics = bam.getBaseName() + "-dedup.stats"

    """
    gatk MarkDuplicatesSpark --java-options "-Xmx${task.memory.giga}g" \
        -I $bam \
        -O $bam_md \
        -M $metrics \
        -OBI --tmp-dir . 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
