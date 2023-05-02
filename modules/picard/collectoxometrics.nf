process OXO_METRICS {

    tag "${meta.patient_id}|${meta.sample_id}"

    label 'medium_serial'

    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Processing/Picard_Metrics", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path(targets)
    tuple path(fasta),path(fai),path(dict)
    tuple path(dbsnp),path(dbsnp_tbi)

    output:
    path(outfile)
    path("versions.yml"), emit: versions

    script:
    outfile = "${meta.patient_id}_${meta.sample_id}.OxoG_metrics.txt"

    """

    picard -Xmx${task.memory.toGiga()-1}G CollectOxoGMetrics \
        INPUT=${bam} \
        OUTPUT=${outfile} \
        DB_SNP=${dbsnp} \
        INTERVALS=${targets} \
        REFERENCE_SEQUENCE=${fasta} \
        TMP_DIR=tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS

    """
}

