process CNVKIT_BATCH {

    label 'medium_parallel'

    tag "${meta.patient_id}|${meta.sample_id}|${cnn}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVkit", mode: 'copy'
    
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    tuple val(meta),path(bam),path(bai)
    path(cnn)

    output:
    tuple val(meta),path(results), emit: results
    tuple val(meta),path(cns), emit: cns
    path("versions.yml"), emit: versions

    script:
    results = "cnvkit_${cnn}_${meta.sample_id}"
    cns = results + "/" + bam.getBaseName() + ".call.cns"

    def options = ""
    def pre = ""
    def cnn_clean = ""

    if ( cnn.getName().contains(".gz") ) {
        cnn_clean = cnn.getBaseName()
        pre = "gunzip -c $cnn > $cnn_clean"
    } else {
        cnn_clean = cnn
    }

    if (meta.status == 1) {
        options += " --segment-method ${params.cnvkit_mode_tumor}"
    } else {
        options += " --segment-method ${params.cnvkit_mode}"
    }

    """
    $pre

    cnvkit.py batch $bam -r $cnn_clean $options -d $results -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$( cnvkit.py version )
    END_VERSIONS
    """

}
