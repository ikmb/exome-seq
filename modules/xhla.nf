process XHLA {

	container "humanlongevity/hla:latest"
	
	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/xHLA", mode: 'copy'

	input:
	tuple val(meta),path(bam),path(bai)

	output:
	tuple val(meta),path("results/*.json"), emit: results

	script:
	def sample_id = "${meta.patient_id}_${meta.sample_id}"

	script:

	"""
		run.py --sample_id $sample_id --input_bam $bam --output_path results --delete
	"""
}
