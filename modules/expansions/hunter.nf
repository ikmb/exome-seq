process EXPANSION_HUNTER {

	container 'quay.io/biocontainers/expansionhunter:4.0.2--he785bd8_0'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/ExpansionHunter", mode: 'copy'

        input:
        tuple val(meta), path(bam),path(bai)
        path(catalog)

        output:
        tuple val(meta),path(expansion_report)
        path(expansion_vcf)

        script:
        expansion_report = "${meta.patient_id}_${meta.sample_id}.expansion_report.json"
        expansion_vcf = "${meta.patient_id}_${meta.sample_id}.expansion_report.vcf"
        prefix = "${meta.patient_id}_${meta.sample_id}.expansion_report"

        """
                ExpansionHunter --reads $bam --reference $params.fasta --variant-catalog $catalog --output-prefix $prefix
        """

}

