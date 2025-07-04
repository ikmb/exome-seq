process CNVKIT_BATCH_PAIRED {

    label 'medium_parallel'

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVkit", mode: 'copy'
    
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    tuple val(meta),path(bam),path(bai),path(tbam),path(tbai)
    path(bed)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(results), emit: results
    tuple val(meta),path("${results}/*dedup.call.cns"), optional: true, emit: cns
    path("versions.yml"), emit: versions

    script:
    results = "cnvkit_${meta.sample_id}"

    """
    cnvkit.py batch $tbam  -n $bam -t $bed -d $results -p ${task.cpus} --segment-method ${params.cnvkit_mode_tumor}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$( cnvkit.py version )
    END_VERSIONS
    """

}
