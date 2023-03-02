process GATK_GET_PILEUP_SUMMARIES {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai)
	path(intervals)

	output:
	tuple val(meta),path(stable), emit: table

	script:
	stable = bam.getBaseName() + "_pileup-summaries.table"

	"""
		gatk GetPileupSummaries \
			-I $bam \
			-V $params.gnomad_af_vcf \
			-L $intervals \
			-O $stable
	"""
	
}
