process GATK_FILTER_MUTECT_CALLS {

	label 'gatk'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi),path(vcf_stats),path(read_orientation_model),path(contamination_table)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf_filtered),path(tbi_filtered), emit: vcf

	script:

	vcf_filtered = vcf.getBaseName() + ".filtered.vcf.gz"
	tbi_filtered = vcf_filtered + ".tbi"

	"""
		gatk FilterMutectCalls \
			-V $vcf \
			-R $fasta \
			-O $vcf_filtered \
			--ob-priors $read_orientation_model \
			--contamination-table $contamination_table \
			-OVI true
			
	"""
}
