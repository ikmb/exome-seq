process GATK_CALCULATE_CONTAMINATION {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

        container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'
    
        label 'medium_serial'

	input:
	tuple val(meta),path(stable)

	output:
	tuple val(meta),path(ctable), emit: table
	path("versions.yml"), emit: versions

	script:
	ctable = stable.getBaseName() + "-contamination.table"

    """
    gatk CalculateContamination \
        -I $stable \
        -O $ctable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
