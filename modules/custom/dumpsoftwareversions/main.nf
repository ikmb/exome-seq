process CUSTOM_DUMPSOFTWAREVERSIONS {

    label 'short_serial'

    container 'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    script:
    template 'dumpsoftwareversions.py'

}
