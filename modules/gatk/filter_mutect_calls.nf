process GATK_FILTER_MUTECT_CALLS {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

	label 'medium_serial'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi),path(vcf_stats),path(read_orientation_model),path(contamination_table)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf_filtered),path(tbi_filtered), emit: vcf
	path("versions.yml"), emit: versions

	script:

	vcf_filtered = vcf.getSimpleName() + "-filtered.vcf.gz"
	tbi_filtered = vcf_filtered + ".tbi"

    """
    gatk FilterMutectCalls \
        -V $vcf \
        -R $fasta \
        -O $vcf_filtered \
        --ob-priors $read_orientation_model \
        --contamination-table $contamination_table \
        -OVI true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
			
    """
}
