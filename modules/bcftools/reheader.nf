process BCFTOOLS_REHEADER {

    tag "${meta.assembly}|${meta.id}"

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    //publishDir "${params.outdir}/${assembly}", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(tbi)
    tuple val(meta_f),path(fai)

    output:
    tuple val(meta),path(vcf_clean),path(tbi_clean), emit: vcf
    path("versions.yml"), emit: versions

    script:
    assembly = meta.assembly
    vcf_clean = vcf.getSimpleName() + ".c.vcf.gz"
    tbi_clean = vcf_clean + ".tbi"

    """
    awk -F"\t" '{ print \$1,1, \$2}' OFS="\t" $fai > regions.txt

    bcftools reheader --fai $fai -o tmp.vcf.gz $vcf
    bcftools index tmp.vcf.gz
    bcftools view -R regions.txt -O z -o $vcf_clean tmp.vcf.gz
    
    bcftools index -t $vcf_clean
    rm tmp.vcf.gz*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}

