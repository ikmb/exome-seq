process BCFTOOLS_SORT {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    input:
    tuple val(meta),path(vcf),path(tbi)

    output:
    tuple val(meta),path(vcf_sorted),path(tbi_sorted), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf_sorted = vcf.getSimpleName() + "_sorted.vcf.gz"
    tbi_sorted = vcf_sorted + ".tbi"

    """    
    bcftools sort -o $vcf_sorted $vcf
    bcftools index -t $vcf_sorted

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}

