process BWA_MEM {

	tag "${meta.patient_id}|${meta.sample_id}"

	//scratch true	

	input:
	tuple val(meta), path(left),path(right)
	val(bwa_index)

	output:
	tuple val(meta), path(bam), emit: bam
	val(sample), emit: sample_name
        val(meta), emit: meta_data
    
	script:
	bam = "${meta.sample_id}_${meta.library_id}_${meta.readgroup_id}.bwa-aligned.fm.bam"
	sample = "${meta.patient_id}_${meta.sample_id}"

	def aligner = "bwa"
	def options = ""
	if (params.bwa2) {
		aligner = "bwa-mem2"
		options = "-K 1000000"
	}
	"""
		$aligner mem $options -H ${params.dict} -M -R "@RG\\tID:${meta.readgroup_id}\\tPL:ILLUMINA\\tPU:${meta.platform_unit}\\tSM:${meta.patient_id}_${meta.sample_id}\\tLB:${meta.library_id}\\tDS:${bwa_index}\\tCN:${meta.center}" \
			-t ${task.cpus} ${bwa_index} $left $right \
			| samtools fixmate -@ ${task.cpus} -m - - \
			| samtools sort -@ ${task.cpus} -m 4G -O bam -o $bam - 
	"""	
}
