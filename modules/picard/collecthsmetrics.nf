process HYBRID_CAPTURE_METRICS {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Processing/Picard_Metrics", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(targets)
    path(baits)
    tuple path(fasta),path(fai),path(dict)

    output:
    path(outfile)
    path(outfile_per_target)
    path("versions.yml"), emit: versions

    script:
    outfile = "${meta.patient_id}_${meta.sample_id}-hybrid_selection_metrics.txt"
    outfile_per_target = "${meta.patient_id}_${meta.sample_id}-hybrid_selection_per_target_metrics.txt"

    """
    picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
        INPUT=${bam} \
        OUTPUT=${outfile} \
        PER_TARGET_COVERAGE=${outfile_per_target} \
        TARGET_INTERVALS=${targets} \
        BAIT_INTERVALS=${baits} \
        REFERENCE_SEQUENCE=${fasta} \
        MINIMUM_MAPPING_QUALITY=$params.min_mapq \
        TMP_DIR=tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0 
    END_VERSIONS
    """
}

