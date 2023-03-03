process GATK_FILTERVARIANTTRANCHES {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'
	
	input:
	tuple val(meta),path(vcf),path(tbi)
	path(snps)
	path(snps_tbi)
	path(indels)
	path(indels_tbi)

	output:
	tuple val(meta),path(vcf_filtered),path(vcf_filtered_tbi), emit: vcf

	script:
	
	vcf_filtered = vcf.getSimpleName() + "-filtered_tranches.vcf.gz"
	vcf_filtered_tbi = vcf_filtered + ".tbi"

	"""
		gatk FilterVariantTranches \
			-V $vcf \
			--resource ${snps.join(' --resource ')} \
			--resource ${indels.join(' --resource ')} \
			--info-key CNN_1D \
			--invalidate-previous-filters \
			--snp-tranche 99.9 --indel-tranche 99.9 \
			-O $vcf_filtered

		gatk IndexFeatureFile -I $vcf_filtered
	"""

}
