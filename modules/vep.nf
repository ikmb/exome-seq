process VEP {

	container 'quay.io/biocontainers/ensembl-vep:99.2--pl526hecc5488_0'

        tag "${meta.patient_id}|${meta.sample_id}"
        
        publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(vcf_index)
	tuple path(fasta),path(fai),path(dict)

        output:
        path(vcf_vep)
        path(vcf_alissa)
        path('*.html')

        script:
        vcf_vep = vcf.getSimpleName() + "-vep.vcf"
	vcf_alissa = vcf.getSimpleName() + "-vep2alissa.vcf"

        """
                vep --offline \
                        --cache \
                        --dir ${params.vep_cache_dir} \
                        --species homo_sapiens \
                        --assembly GRCh38 \
                        -i $vcf \
                        --format vcf \
                        -o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
                        --plugin dbNSFP,${params.dbnsfp_db},${params.dbnsfp_fields} \
                        --plugin dbscSNV,${params.dbscsnv_db} \
                        --plugin CADD,${params.cadd_snps},${params.cadd_indels} \
                        --plugin ExACpLI \
                        --plugin UTRannotator \
                        --plugin Mastermind,${params.vep_mastermind} \
                        --plugin SpliceAI,${params.spliceai_fields} \
			--af_gnomad \
                        --fasta $fasta \
                        --fork ${task.cpus} \
                        --vcf \
                        --per_gene \
                        --sift p \
                        --polyphen p \
                        --check_existing \
                        --canonical

                sed -i.bak 's/CADD_PHRED/CADD_phred/g' $vcf_vep
                vep2alissa.pl --infile $vcf_vep > $vcf_alissa
        """

}

