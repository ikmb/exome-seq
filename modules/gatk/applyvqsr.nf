process GATK_APPLYVQSR {

    tag "${meta.patient_id}|${meta.sample_id}"
    
    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    //publishDir "${params.outdir}/GATK", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)
    tuple path(recal),path(recal_idx)
    path(tranches)
    val(modus)

    output:
    tuple val(meta),path(vcf_recal),path(vcf_recal_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_recal = vcf.getSimpleName() + "-" + modus + "-recal.vcf.gz"
    vcf_recal_tbi = vcf_recal + ".tbi"

    """
    gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
        -V $vcf \
        --recal-file $recal \
        --tranches-file $tranches \
        --truth-sensitivity-filter-level 99.7 \
        --create-output-variant-index true \
        -mode $modus \
        -O $vcf_recal

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
