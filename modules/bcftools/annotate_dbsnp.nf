process BCFTOOLS_ANNOTATE_DBSNP {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    tag "${meta.patient_id}|${meta.sample_id}|${meta.variantcaller}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)
    tuple path(dbsnp),path(dbsnp_tbi)

    output:
    tuple val(meta),path(vcf_annotated), path(vcf_annotated_index), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_annotated = vcf.getSimpleName() + "-rsids.vcf.gz"
    vcf_annotated_index = vcf_annotated + ".tbi"

    """
    bcftools annotate -c ID -a $dbsnp -O z -o $vcf_annotated $vcf
    tabix $vcf_annotated

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """
}

