process MULTIQC {

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    publishDir "${params.outdir}/Summary/${cname}", mode: 'copy'

    input:
    val(cname)
    path('*')

    output:
    path("${cname}_multiqc.html"), emit: report
    path("versions.yml"), emit: versions

    script:

    """
    cp $params.logo .
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n ${cname}_multiqc *

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS

    """
}

process MULTIQC_PANEL {

    container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

    publishDir "${params.outdir}/Summary/Panel", mode: "copy"

    input:
    tuple val(panel_name),path('*')

    output:
    path("${panel_name}_multiqc.html")

    script:

    """
    cp $params.logo .
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n ${panel_name}_multiqc *

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

