process GATK_CNNSCOREVARIANTS {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'

	container = "broadinstitute/gatk:4.2.6.1"
	
	label 'medium_serial'

	input:
	tuple val(meta),path(vcf),path(tbi),path(bam),path(bai)
	path(intervals)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(cnn_vcf),path(cnn_vcf_tbi), emit: vcf
	path("versions.yml"), emit: versions
	
	script:
	
	cnn_vcf = vcf.getSimpleName() + "-cnn.vcf.gz"
	cnn_vcf_tbi = cnn_vcf + ".tbi"

	"""
    gatk CNNScoreVariants \
        --variant $vcf \
        --output $cnn_vcf \
        --reference $fasta \
        --intervals $intervals \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
	"""
}
