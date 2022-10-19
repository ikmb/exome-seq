process GATK_MERGEVCFS {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/VQSR", mode: 'copy'	

	label 'gatk'

	input:
	tuple val(meta),path(vcf_snp),path(vcf_indel)
	tuple val(tmeta),path(tbi_snp),path(tbi_indel)

	output:
	tuple val(meta),path(merged_vcf),path(merged_vcf_tbi), emit: vcf

	script:
	merged_vcf = "gatk-" + params.run_name + ".merged.vcf.gz"
	merged_vcf_tbi = merged_vcf + ".tbi"

	"""
		gatk --java-options "-Xmx4g" MergeVcfs \
			--INPUT $vcf_snp --INPUT $vcf_indel \
			--OUTPUT $merged_vcf
		gatk IndexFeatureFile -I $merged_vcf
	"""
}
