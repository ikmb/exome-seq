process vep {

	label 'vep'

        publishDir "${params.outdir}/${indivID}/${sampleID}/VEP/${cname}", mode: 'copy'

        input:
        tuple val(cname),val(indivID),val(sampleID),path(vcf),path(vcf_index)

        output:
        path(vcf_vep)
        path(vcf_alissa)
        path('*.html')

        script:
        vcf_vep = vcf.getBaseName() + ".vep.vcf"
        vcf_alissa = vcf.getBaseName() + ".vep2alissa.vcf"

        """
                vep --offline \
                        --cache \
                        --dir ${params.vep_cache_dir} \
                        --species homo_sapiens \
                        --assembly $params.assembly \
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
                        --fasta $params.fasta \
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
