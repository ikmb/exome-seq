process GLNEXUS {

    container 'quay.io/mlin/glnexus:v1.3.1'

    label 'medium_serial'

    publishDir "${params.outdir}/MERGED_CALLSET/GLNEXUS_DEEPVARIANT", mode: 'copy'

    input:
    path(gvcfs)
    path(bed)

    output:
    tuple path(merged_vcf),path(merged_vcf_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    merged_vcf = "deepvariant-joint_calling-" + params.run_name + ".vcf.gz"
    merged_vcf_tbi = merged_vcf + ".tbi"
    """
    /usr/local/bin/glnexus_cli \
        --config ${params.glnexus_config} \
        --bed $bed \
        $gvcfs | bcftools view - | bgzip -c > $merged_vcf

    tabix $merged_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        glnexus: 1.4.1
    END_VERSIONS

    """
}

