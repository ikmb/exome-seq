process GATK_LEARN_READ_ORIENTATION_MODEL {

	label 'gatk'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

	input:
	tuple val(meta),path(f12r)

	output:
	tuple val(meta),path(model), emit: model

	script:
	model = f12r.getBaseName() + "_read-orientation-model.tar.gz"

	"""
		gatk LearnReadOrientationModel -I $f12r -O $model
	"""

}
