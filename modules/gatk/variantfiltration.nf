process GATK_VARIANTFILTRATION {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'

        container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'short_serial'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple path(vcf_filtered),path(vcf_filtered_tbi), emit: vcf
	path("versions.yml"), emit: versions

	script:
	vcf_filtered = vcf.getBaseName() + "-filtered.vcf.gz"
	vcf_filtered_tbi = vcf_filtered +".tbi"
	def options = ""
	options = "--filter-expression \"${params.gatk_hard_filter}\" --filter-name HardFilter"

    """
    gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
        -V $vcf \
        -O $vcf_filtered -OVI true \
        $options

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
