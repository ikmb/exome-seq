process GATK_SELECTVARIANTS {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'short_serial'

    input:
    tuple val(m_f),path(vcf),path(vcf_index)
    val(meta)

    output:
    tuple val(meta),path(vcf_sample),path(vcf_sample_index), emit: vcf
    path("versions.yml"), emit: versions

    script:
    def prefix = meta.sample_id
    vcf_sample = prefix + "-" + m_f.variantcaller + "-split.vcf.gz"
    vcf_sample_index = vcf_sample + ".tbi"

    """
    gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $prefix -O $vcf_sample -OVI

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

    """

}

