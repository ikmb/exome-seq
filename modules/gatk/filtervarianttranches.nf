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
			--info-key CNN_2D \
			--snp-tranche 99.9 --snp-tranche 99.95 \
			--indel-tranche 99.0 --indel-tranche 99.4 \
			--invalidate-previous-filters \
			-O $vcf_filtered

		gatk IndexFeatureFile -I $vcf_filtered
	"""

}
