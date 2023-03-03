process GATK_MARK_DUPLICATES {

	label 'gatk'

	tag "${meta.patient_id}|${meta.sample_id}"

	input:
	tuple val(meta),path(bam),path(bai)

	output:
	tuple val(meta),path(bam_md),path(bai_md), emit: bam
	tuple val(meta),path(stats), emit: stats

	script:

	bam_md = bam.getBaseName() + "-dedup.bam"
	bai_md = bam.getBaseName() + "-dedup.bam.bai"
	metrics = bam.getBaseName() + "-dedup.stats"

	"""
		gatk MarkDuplicatesSpark --java-options "-Xmx${task.memory.giga}g" \
			-I $bam \
			-O $bam_md \
			-M $metrics \
			-OBI --tmp-dir . 
	"""

}
