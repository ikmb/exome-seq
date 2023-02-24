process GATK_MUTECT2 {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf),path(tbi), emit: vcf

	script:
	vcf = bam.getBaseName() + ".mutect2.vcf.gz"
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
		-I $bam \
		-O $vcf \
		-L $intervals \
		-OVI \
	"""	

}
