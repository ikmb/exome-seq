process FASTQ_TO_BAM {

	container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
	
	tag "${meta.patient_id}|${meta.sample_id}"

	input:
	tuple val(meta),path(R1),path(R2)

	output:
	tuple val(meta),path(bam), emit: bam
	path("versions.yml"), emit: versions

	script:

	bam = meta.patient_id + "_" + meta.sample_id + "_" + meta.library_id + "-u.bam"

	"""
	picard FastqToSam F1=$R1 F2=$R2 O=$bam SM=${meta.sample_id} RG=${meta.readgroup_id} LB=${meta.library_id} R=${params.fasta} PU=${meta.platform_unit}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS
    """

}

	
