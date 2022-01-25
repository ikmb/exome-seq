process align {

	//scratch true	

	input:
	tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(run_date), val(center),path(left),path(right)

	output:
	tuple val(indivID), val(sampleID), path(bam)
    
	script:
	bam = sampleID + "_" + libraryID + "_" + rgID + ".aligned.fm.bam"	

	def aligner = "bwa"
	def options = ""
	if (params.bwa2) {
		aligner = "bwa-mem2"
		options = "-K 1000000"
	}
	"""
		$aligner mem $options -H ${params.dict} -M -R "@RG\\tID:${rgID}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tSM:${indivID}_${sampleID}\\tLB:${libraryID}\\tDS:${params.fasta}\\tCN:${center}" \
			-t ${task.cpus} ${params.bwa_index} $left $right \
			| samtools fixmate -@ ${task.cpus} -m - - \
			| samtools sort -@ ${task.cpus} -m 3G -O bam -o $bam - 
	"""	
}
