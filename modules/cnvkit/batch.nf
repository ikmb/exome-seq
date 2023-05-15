process CNVKIT_BATCH {

    label 'medium_parallel'

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVkit", mode: 'copy'
    
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    tuple val(meta),path(bam),path(bai)
    path(cnn)
    path(targets)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(results), emit: results
    tuple val(meta),path(cns), emit: cns
    path("versions.yml"), emit: versions

    script:
    results = "cnvkit_${meta.sample_id}"
    cns = results + "/" + bam.getBaseName() + ".call.cns"

    def options = ""
    def pre = ""
    def cnn_clean = ""

    if (cnn) {
        if ( cnn.getName().contains(".gz") ) {
           cnn_clean = cnn.getBaseName()
           pre = "gunzip -c $cnn > $cnn_clean"
        }
        options += " -r $cnn_clean"
    } else {
        options += " -f $fasta "
        if (targets) {
            options += " -t $targets"
        }
    }

    if (meta.status == 1) {
        options += " --segment-method ${params.cnvkit_mode_tumor}"
    } else {
        options += " --segment-method ${params.cnvkit_mode}"
    }

    """
    $pre

    cnvkit.py batch $bam $options -d $results -p ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$( cnvkit.py version )
    END_VERSIONS
    """

}
