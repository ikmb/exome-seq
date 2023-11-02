process GATK_GATHERPILEUPSUMMARIES {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    input:
    tuple val(meta),path(tables)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(merged_table), emit: table
    path("versions.yml"), emit: versions

    script:
    merged_table = meta.sample_id + "-pileup_summaries.table"

    """
    gatk GatherPileupSummaries \
        -SD $dict \
        -I ${tables.join(' -I ')} \
        -O $merged_table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
}
