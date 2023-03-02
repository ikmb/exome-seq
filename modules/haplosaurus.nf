process HAPLOSAURUS {

	container 'quay.io/biocontainers/ensembl-vep:99.2--pl526hecc5488_0'

        tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(vcf_index)
	tuple path(fasta),path(fai),path(dict)

        output:
        path(haplo)

        script:
        haplo = vcf.getBaseName() + ".vep_haplo.txt"

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

