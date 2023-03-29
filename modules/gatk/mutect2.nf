process GATK_MUTECT2 {

	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf),path(tbi),path(stats), emit: vcf
	tuple val(meta),path(f1r2), emit: f1r2
	path("versions.yml"), emit: versions

	script:
	vcf = bam.getBaseName() + "-mutect2.vcf.gz"
	tbi = vcf + ".tbi"
	stats = vcf + ".stats"
	f1r2 =  bam.getBaseName() + "-f1r2.tar.gz"

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
        --f1r2-tar-gz $f1r2 \
        -OVI \
        $options

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """	

}
