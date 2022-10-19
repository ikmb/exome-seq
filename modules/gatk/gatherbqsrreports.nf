process GATK_GATHERBQSRREPORTS {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	input:
	tuple val(meta),path(reports)

	output:
	tuple val(meta),path(merged_report), emit: report

	script:
	merged_report = meta.patient_id + "_" + meta.sample_id + "_recal.txt"

	"""
		gatk GatherBQSRReports --input ${reports.join(' --input ')} --output $merged_report
	"""

}
