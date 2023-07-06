process GATK_LEARN_READ_ORIENTATION_MODEL {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'short_serial'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

    input:
    tuple val(meta),path(f12rs)

    output:
    tuple val(meta),path(model), emit: model
    path("versions.yml"), emit: versions

    script:
    model = meta.sample_id + "-read_orientation_model.tar.gz"

    """
    gatk LearnReadOrientationModel -I ${f12rs.join(' -I ')} -O $model

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
