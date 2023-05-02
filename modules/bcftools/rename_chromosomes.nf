process BCFTOOLS_RENAME_CHROMOSOMES {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    publishDir "${params.outdir}/${assembly}/", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)
    path(lookup)

    output:
    tuple val(meta),path(vcf_r),path(tbi_r), emit: vcf
    path("versions.yml"), emit: versions

    script:
    assembly = meta.assembly
    vcf_r = vcf.getSimpleName() + "_final.vcf.gz"
    tbi_r = vcf_r + ".tbi"

    """
    echo "##reference=${params.assembly}" > header.txt
    bcftools annotate --rename-chrs ${lookup} -h header.txt -O z -o $vcf_r $vcf
    tabix $vcf_r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}


