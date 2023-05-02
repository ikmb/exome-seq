process GATK_MERGEVCFS {

	tag "ALL"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/VQSR", mode: 'copy'	

        container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'medium_serial'

	input:
	tuple val(meta),path(vcf_snp),path(vcf_indel)
	tuple val(tmeta),path(tbi_snp),path(tbi_indel)

	output:
	tuple val(meta),path(merged_vcf),path(merged_vcf_tbi), emit: vcf
	path("versions.yml"), emit: versions

	script:
	merged_vcf = "gatk-" + params.run_name + "-merged.vcf.gz"
	merged_vcf_tbi = merged_vcf + ".tbi"

    """
    gatk --java-options "-Xmx4g" MergeVcfs \
        --INPUT $vcf_snp --INPUT $vcf_indel \
        --OUTPUT $merged_vcf
    gatk IndexFeatureFile -I $merged_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
	"""
}
