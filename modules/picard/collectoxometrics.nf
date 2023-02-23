process OXO_METRICS {

        tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/picard:2.20.4--1'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Processing/Picard_Metrics", mode: 'copy'

        input:
        tuple val(meta), path(bam), path(bai)
        path(targets)
	tuple path(fasta),path(fai),path(dict)
	tuple path(dbsnp),path(dbsnp_tbi)

        output:
        path(outfile)

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
        """
}

