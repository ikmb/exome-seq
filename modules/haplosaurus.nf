process HAPLOSAURUS {

    label 'medium_parallel'

    container 'quay.io/biocontainers/ensembl-vep:109.3--pl5321h4a94de4_0'

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(vcf_index)
    tuple path(fasta),path(fai),path(dict)

    output:
    path(haplo)
    path("versions.yml"), emit: versions

    script:
    haplo = vcf.getSimpleName() + "-vep_haplo.txt"

    """
    haplo --offline \
        --cache \
        --dir ${params.vep_cache_dir} \
        --species homo_sapiens \
        --assembly GRCh38 \
        -i $vcf \
        --format vcf \
        -o $haplo  \
        --fasta $fasta \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplosaurus: \$( echo \$(haplo --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS

    """
}

