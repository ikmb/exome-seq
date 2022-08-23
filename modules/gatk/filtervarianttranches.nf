process GATK_FILTERVARIANTTRANCHES {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'
	
	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_filtered),path(vcf_filtered_tbi), emit: vcf

	script:
	
	vcf_filtered = vcf.getBaseName() + ".filtered_tranches.vcf.gz"
	vcf_filtered_tbi = vcf_filtered + ".tbi"

	"""
		gatk FilterVariantTranches \
			-V $vcf \
			--resource $params.hapmap \
			--resource $params.mills \
			--resource $params.omni \
			--resource $params.g1k \
			--resource $params.dbsnp \
			--info-key CNN_2D \
			--invalidate-previous-filters \
			--snp-tranche 99.9 --indel-tranche 99.9 \
			-O $vcf_filtered

		gatk IndexFeatureFile -I $vcf_filtered
	"""

}
