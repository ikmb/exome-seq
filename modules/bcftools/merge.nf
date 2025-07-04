process BCFTOOLS_MERGE {

    container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

    publishDir "${params.outdir}/MERGED_CALLSET/BCFTOOLS", mode: 'copy'

    input:
    tuple val(meta),path(vcfs),path(tbis)

    output:
    tuple val(meta),path(merged_vcf),path(merged_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    merged_vcf = meta.variantcaller + ".flat_merged." + params.run_name + ".vcf.gz"
    merged_tbi = merged_vcf + ".tbi"
    
    """
    bcftools merge --threads ${task.cpus} -o $merged_vcf -O z *.vcf.gz
    bcftools index -t $merged_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version | head -n1 | sed -e "s/bcftools //g" )
    END_VERSIONS
    """

}

