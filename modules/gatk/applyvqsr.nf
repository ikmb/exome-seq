process GATK_APPLYVQSR {

	label 'gatk'

	//publishDir "${params.outdir}/GATK", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi)
	tuple path(recal),path(recal_idx)
	path(tranches)
	val(modus)

	output:
	tuple val(meta),path(vcf_recal),path(vcf_recal_tbi), emit: vcf

	script:
	vcf_recal = vcf.getSimpleName() + "-" + modus + "-recal.vcf.gz"
	vcf_recal_tbi = vcf_recal + ".tbi"

	"""
		gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
			-V $vcf \
			--recal-file $recal \
			--tranches-file $tranches \
			--truth-sensitivity-filter-level 99.7 \
			--create-output-variant-index true \
			-mode $modus \
			-O $vcf_recal
	"""

}
