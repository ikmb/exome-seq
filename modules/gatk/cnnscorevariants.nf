process GATK_CNNSCOREVARIANTS {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(vcf),path(tbi),path(bam),path(tbi)
	path(intervals)

	output:
	tuple val(meta),path(cnn_vcf),path(cnn_vcf_tbi), emit: vcf
	
	script:
	
	cnn_vcf = vcf.getBaseName() + ".cnn.vcf.gz"
	cnn_vcf_tbi = cnn_vcf + ".tbi"

	"""
		gatk CNNScoreVariants \
			--variant $vcf \
			--input $bam \
			--output $cnn_vcf \
			--reference $params.fasta \
			--intervals $intervals \
			
	"""
	
}
