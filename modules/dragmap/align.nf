process DRAGMAP_ALIGN {

	container 'quay.io/biocontainers/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:5ebebbc128cd624282eaa37d2c7fe01505a91a69-0'

	input:
	tuple val(meta),path(R1),path(R2)
	val(dragen_ref)

	output:
	tuple val(meta),path(bam), emit: bam
	val(sample), emit: sample_name
	path(log_file), emit: log

	script:	
	bam = "${meta.sample_id}_${meta.library_id}_${meta.readgroup_id}-dragmap_aligned-fm.bam"
        sample = "${meta.patient_id}_${meta.sample_id}"
	log_file = sample + "-dragmap.txt"
	"""
		dragen-os \
			-r $dragen_ref \\
			-1 $R1 \\
			-2 $R2 \\
			--RGID ${meta.readgroup_id} \\
			--RGSM $sample \\
			--num-threads ${task.cpus} 2> $log_file | samtools fixmate -m --threads ${task.cpus} - - | samtools sort --threads ${task.cpus} -m 2G -O bam -o $bam -

	"""
}