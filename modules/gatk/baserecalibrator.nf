process GATK_BASERECALIBRATOR {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai),path(intervals)
	tuple path(fasta),path(fai),path(dict)
	path(snps)
	path(snps_tbi)
	path(indels)
	path(indels_tbi)

	output:
	tuple val(meta),path(report), emit: report

	script:
	report = bam.getBaseName() + "-" + intervals.getBaseName() + "_recal.txt"	

	"""
		gatk BaseRecalibrator -R ${fasta}  -I $bam -O $report \
			-L $intervals --use-original-qualities --known-sites ${snps.join(' --known-sites ')} \
			--known-sites ${indels.join(' --known-sites ')} \
			-ip $params.interval_padding
	"""
}
