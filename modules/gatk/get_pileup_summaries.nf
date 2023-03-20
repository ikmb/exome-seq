process GATK_GET_PILEUP_SUMMARIES {

	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(stable), emit: table

	script:
	stable = bam.getBaseName() + "-pileup_summaries.table"

	"""
		gatk GetPileupSummaries \
			-I $bam \
			-V $params.gnomad_af_vcf \
			-L $intervals \
			-R $fasta \
			-O $stable
	"""
	
}
