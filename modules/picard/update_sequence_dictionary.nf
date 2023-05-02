process PICARD_UPDATE_SEQUENCE_DICTIONARY {

    tag "${meta.id}"

    container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'

    publishDir "${params.outdir}/${assembly}", mode: 'copy'

    label 'short_serial'

    input:
    tuple val(meta),path(vcf),path(tbi)
    tuple val(meta_d),path(dict)

    output:
    tuple val(meta),path(vcf_update),path(tbi_update), emit: vcf
    path("versions.yml"), emit: versions

    script:
    assembly = meta.assembly
    vcf_update = meta.file_name
    tbi_update = vcf_update + ".tbi"

    """
    picard UpdateVcfSequenceDictionary -I $vcf -O $vcf_update -SD $dict --CREATE_INDEX true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS

    """

}

