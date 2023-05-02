process VEP_VEP {

    label 'medium_parallel'

    container 'quay.io/biocontainers/ensembl-vep:109.3--pl5321h4a94de4_0'

    tag "${meta.patient_id}|${meta.sample_id}"
        
    publishDir "${params.outdir}/VEP/${meta.variantcaller}", mode: 'copy'

    input:
    tuple val(meta),path(vcf),path(vcf_index)
    tuple path(fasta),path(fai),path(dict)

    output:
    path(vcf_vep)
    path(vcf_alissa)
    path('*.html')
    path("versions.yml"), emit: versions

    script:
    vcf_vep = vcf.getSimpleName() + "-vep.vcf"
    vcf_alissa = vcf.getSimpleName() + "-vep2alissa.vcf"

    def options = ""
    if (params.dbnsfp_db) {
        options.concat(" --plugin dbNSFP,${params.dbnsfp_db},${params.dbnsfp_fields}")
    }
    if (params.dbscsnv_db) {
        options.concat("--plugin dbscSNV,${params.dbscsnv_db}")
    }
    if (params.cadd_snps) {
        options.concat("--plugin CADD,${params.cadd_snps},${params.cadd_indels}") 
    }
    if (params.vep_mastermind) {
        options.concat("--plugin Mastermind,${params.vep_mastermind}")
    }
    if (params.spliceai_fields) {
        options.concat("--plugin SpliceAI,${params.spliceai_fields}")
    }

    """
    vep --offline \
        --cache \
        --dir ${params.vep_cache_dir} \
        --species homo_sapiens \
        --assembly GRCh38 \
        -i $vcf \
        --format vcf \
        -o $vcf_vep --dir_plugins ${params.vep_plugin_dir} \
        $options \
        --plugin ExACpLI \
        --plugin UTRannotator \
        --af_gnomade \
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | grep ensembl-vep | sed 's/  ensembl-vep          ://')
    END_VERSIONS
    """

}

