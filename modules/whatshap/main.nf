process WHATSHAP {

	container 'quay.io/biocontainers/whatshap:1.1--py36hae55d0a_1'

	tag "${meta.patient_id}|${meta.sample_id}"
		
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi)
	path(bams)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(phased_vcf),path(phased_tbi), emit: vcf

	script:
	phased_vcf = vcf.getSimpleName() + "-phased.vcf.gz"
	phased_tbi = phased_vcf + ".tbi"

	"""
		whatshap phase -o $phased_vcf --tag=PS --reference $fasta $vcf *.cram
		tabix $phased_vcf
	"""
		
}

process WHATSHAP_SINGLE {

	tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/whatshap:1.1--py36hae55d0a_1'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi),path(bam),path(bai)
	tuple path(fasta),path(fai),path(dict)

        output:
        tuple val(meta),path(phased_vcf),path(phased_tbi), emit: vcf

        script:
        phased_vcf = vcf.getSimpleName() + "-phased.vcf.gz"
        phased_tbi = phased_vcf + ".tbi"

        """
                whatshap phase -o $phased_vcf --tag=PS --reference $fasta $vcf *.cram
                tabix $phased_vcf
        """

}

