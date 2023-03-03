process GATK_MUTECT2_MULTI {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2", mode: 'copy'

	label 'gatk'

	input:
	val(meta)
	path(bams)
	path(bais)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple path(vcf),path(tbi), emit: vcf

	script:
	vcf = params.run_name + ".mutect2.joint-calling.vcf.gz"
	tbi = vcf + ".tbi"

	def options = ""
	if (params.mutect_normals) {
		options += " --panel-of-normals ${params.mutect_normals}"
	}
	if (params.gnomad_af_vcf) {
		options += " --germline-resource ${params.gnomad_af_vcf}"
	}

	"""
		gatk Mutect2 \
		-R $fasta \
		-I ${bams.join(' -I ')} \
		-O $vcf \
		-L $intervals \
		-OVI \
		$options
	"""	

}
