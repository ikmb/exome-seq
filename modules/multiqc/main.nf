process multiqc {

    publishDir "${params.outdir}/Summary/#{cname}", mode: 'copy'

    input:
    val(cname)
    path('*')

    output:
    file("${cname}_multiqc.html")
    file("*.yaml")

    script:

    """
    cp $params.logo .
    cp $baseDir/conf/multiqc_config.yaml multiqc_config.yaml
    multiqc -n ${cname}_multiqc *

    """
}

process multiqc_panel {

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
        """
}

