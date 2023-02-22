process HAPLOSAURUS {

	container 'quay.io/biocontainers/ensembl-vep:99.2--pl526hecc5488_0'

        publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(vcf_index)
	tuple path(fasta),path(fai),path(dict)

        output:
        path(vcf_vep)
        path('*.html')

        script:
        vcf_vep = vcf.getBaseName() + ".vep_haplo.txt"

        """
                haplo --offline \
                        --cache \
                        --dir ${params.vep_cache_dir} \
                        --species homo_sapiens \
                        --assembly $params.assembly \
                        -i $vcf \
                        --format vcf \
                        -o $vcf_vep
                        --fasta $fasta \
                        --fork ${task.cpus} \
                        --vcf \
                        --per_gene \
                        --sift p \
                        --polyphen p \
                        --check_existing \
                        --canonical

        """


}

