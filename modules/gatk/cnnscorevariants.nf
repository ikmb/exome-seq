process GATK_CNNSCOREVARIANTS {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(vcf),path(tbi),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(cnn_vcf),path(cnn_vcf_tbi), emit: vcf
	
	script:
	
	cnn_vcf = vcf.getSimpleName() + "-cnn.vcf.gz"
	cnn_vcf_tbi = cnn_vcf + ".tbi"

	"""
		gatk CNNScoreVariants \
			--variant $vcf \
			--output $cnn_vcf \
			--reference $fasta \
			--intervals $intervals \
	"""
	
}
