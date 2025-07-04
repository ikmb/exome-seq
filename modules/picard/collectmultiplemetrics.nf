process MULTI_METRICS {

    tag "${meta.patient_id}|${meta.sample_id}"

    label 'medium_serial'

    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Processing/Picard_Metrics", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(baits)
    tuple path(fasta),path(fai),path(dict)
    tuple path(dbsnp),path(dbsnp_tbi)

    output:
    file("${prefix}*")
    path("versions.yml"), emit: versions

    script:
    prefix = "${meta.patient_id}_${meta.sample_id}"

    """
    picard -Xmx5g CollectMultipleMetrics \
        PROGRAM=MeanQualityByCycle \
        PROGRAM=QualityScoreDistribution \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics\
        PROGRAM=CollectSequencingArtifactMetrics \
        PROGRAM=CollectQualityYieldMetrics \
        PROGRAM=CollectGcBiasMetrics \
        PROGRAM=CollectBaseDistributionByCycle \
        INPUT=${bam} \
        REFERENCE_SEQUENCE=${fasta} \
        DB_SNP=${dbsnp} \
        INTERVALS=${baits} \
        ASSUME_SORTED=true \
        QUIET=true \
        OUTPUT=${prefix} \
        TMP_DIR=tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS
    """
}

