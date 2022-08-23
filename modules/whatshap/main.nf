process WHATSHAP {

	tag "${meta.patient_id}|${meta.sample_id}"
		
	label 'whatshap'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi)
	path(bams)

	output:
	tuple val(meta),path(phased_vcf),path(phased_tbi), emit: vcf

	script:
	phased_vcf = vcf.getSimpleName() + "_phased.vcf.gz"
	phased_tbi = phased_vcf + ".tbi"

	"""
		whatshap phase -o $phased_vcf --tag=PS --reference $params.fasta $vcf *.bam
		tabix $phased_vcf
	"""
		
}

process WHATSHAP_SINGLE {

	tag "${meta.patient_id}|${meta.sample_id}"

        label 'whatshap'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi),path(bam),path(bai)

        output:
        tuple val(meta),path(phased_vcf),path(phased_tbi), emit: vcf

        script:
        phased_vcf = vcf.getSimpleName() + "_phased.vcf.gz"
        phased_tbi = phased_vcf + ".tbi"

        """
                whatshap phase -o $phased_vcf --tag=PS --reference $params.fasta $vcf *.bam
                tabix $phased_vcf
        """

}

