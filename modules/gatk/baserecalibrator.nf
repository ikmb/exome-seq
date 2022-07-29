process GATK_BASERECALIBRATOR {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai),path(intervals)

	output:
	tuple val(meta),path(report), emit: report

	script:
	report = bam.getBaseName() + "-" + intervals.getBaseName() + "_recal.txt"	

	"""
		gatk BaseRecalibrator -R ${params.fasta}  -I $bam -O $report \
			-L $intervals --use-original-qualities --known-sites ${params.snps.join(' --known-sites ')} \
			--known-sites ${params.indels.join(' --known-sites ')} \
			-ip $params.interval_padding
	"""
}
