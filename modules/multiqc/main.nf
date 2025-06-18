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
    cat <<-EOF > multiqc_config.yaml
    title: "Exome pipeline report"
    subtitle: "Workflow for clinical SNP calling"

    extra_fn_clean_exts:
      - _R1
      - _R2
      - _duplicate_metrics.txt
      - .pass

    report_comment: >
      This report has been generated automatically by the IKMB Diagnostic Exome pipeline.
      For help interpreting the outputs, please see: https://github.com/ikmb/exome-seq

    report_header_info:
      - Contact E-mail: "b.loescher@ikmb.uni-kiel.de"
      - Application Type: "Whole exome sequencing"

    top_modules:
      - 'general_stats'
      - 'picard'
      - 'picard_hsmetrics'
      - 'bcftools'
      - 'sex_check'
    EOF


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
    cat <<EOF > multiqc_config.yaml
    title: "Exome pipeline report"
    subtitle: "Workflow for clinical SNP calling"

    extra_fn_clean_exts:
      - _R1
      - _R2
      - _duplicate_metrics.txt
      - .pass

    report_comment: >
      This report has been generated automatically by the IKMB Diagnostic Exome pipeline.
      For help interpreting the outputs, please see: https://github.com/ikmb/exome-seq

    report_header_info:
      - Contact E-mail: "b.loescher@ikmb.uni-kiel.de"
      - Application Type: "Whole exome sequencing"

    top_modules:
      - 'general_stats'
      - 'picard'
      - 'picard_hsmetrics'
      - 'bcftools'
      - 'sex_check'
    EOF


    multiqc -n ${panel_name}_multiqc *

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

