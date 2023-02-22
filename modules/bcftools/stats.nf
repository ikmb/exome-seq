process BCFTOOLS_STATS {

	container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

        tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Stats", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        path(vcf_stats), emit: stats

        script:
        vcf_stats = vcf.getBaseName() + ".stats"

        """
                bcftools stats $vcf > $vcf_stats
        """

}

