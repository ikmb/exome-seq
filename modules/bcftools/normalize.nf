process BCFTOOLS_NORMALIZE {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    input:
    tuple val(meta),path(vcf),path(tbi)

    output:
    tuple val(meta),path(vcf_norm),path(tbi_norm), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_norm = vcf.getSimpleName() + "_normalized.vcf.gz"
    tbi_norm = vcf_norm + ".tbi"

    """
    bcftools norm -O z -f $params.fasta -o $vcf_norm $vcf
    bcftools index -t $vcf_norm

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """
    
}

