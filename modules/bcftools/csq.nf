process BCFTOOLS_CSQ {

        tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf_in),path(tbi_in)
	tuple path(fasta),path(fai),path(dict)
	path(gtf)

        output:
        tuple val(meta),path(vcf_fixed),path(tbi_fixed), emit: vcf

        script:
        vcf_fixed = vcf_in.getSimpleName() + "_csq.vcf.gz"
        tbi_fixed = vcf_fixed + ".tbi"

        """
                bcftools csq -f $fasta -g $gtf --phase a $vcf_in -o $vcf_fixed
                bcftools index -t $vcf_fixed
        """

}

