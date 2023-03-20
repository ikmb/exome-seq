process HAPLOSAURUS {

	container 'quay.io/biocontainers/ensembl-vep:109.3--pl5321h4a94de4_0'

        tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(vcf_index)
	tuple path(fasta),path(fai),path(dict)

        output:
        path(haplo)

        script:
        haplo = vcf.getSimpleName() + "-vep_haplo.txt"

        """
                haplo --offline \
                        --cache \
                        --dir ${params.vep_cache_dir} \
                        --species homo_sapiens \
                        --assembly $params.assembly \
                        -i $vcf \
                        --format vcf \
                        -o $haplo  \
                        --fasta $fasta \

        """


}

