process GATK_CALCULATE_CONTAMINATION {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(stable)

	output:
	tuple val(meta),path(ctable), emit: table

	script:
	ctable = stable.getBaseName() + ".contamination.table"

	"""
		gatk CalculateContamination \
			-I $stable \
			-O $ctable
	"""

}
