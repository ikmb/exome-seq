process CNVKIT_BATCH {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVkit", mode: 'copy'
    
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    tuple val(meta),path(bam),path(bai)
    path(cnn)

    output:
    tuple val(meta),path(results), emit: results
    tuple val(meta),path(cns), emit: cns


    script:
    results = "cnvkit_${meta.sample_id}"
    cns = results + "/" + bam.getBaseName() + ".call.cns"

    """
        cnvkit.py batch -r $cnn $bam -d $results -p ${task.cpus} --segment-method ${params.cnvkit_mode}
    """

}