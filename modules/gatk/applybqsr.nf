process GATK_APPLYBQSR {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai),path(recal)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(recal_bam),path(recal_bai), emit: bam

	script:
	recal_bam = bam.getBaseName() + "-recal.bam"
	recal_bai = bam.getBaseName() + "-recal.bai"

	"""
		gatk ApplyBQSR -R $fasta -I $bam -O $recal_bam -L $intervals -bqsr $recal \
                        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                        --add-output-sam-program-record \
                        --create-output-bam-md5 \
                        --use-original-qualities -OBI true -ip $params.interval_padding

	"""
}
