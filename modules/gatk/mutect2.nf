process GATK_MUTECT2 {

	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

        container 'quay.io/biocontainers/gatk4:4.1.8.0--py38h37ae868_0'

        label 'medium_serial'

	input:
	tuple val(meta),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)
	path(mutect_normals)
	path(mutect_normals_tbi)

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
	if (mutect_normals) {
		options += " --panel-of-normals ${mutect_normals}"
	}

	if (params.gnomad_af_vcf) {
		options += " --germline-resource ${params.gnomad_af_vcf}"
	}

    """

    gatk  --java-options "-Xmx${task.memory.giga}g" Mutect2 \
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
