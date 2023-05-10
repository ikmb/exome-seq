process BWA_MEM {

    tag "${meta.patient_id}|${meta.sample_id}"

    label 'medium_parallel'

    input:
    tuple val(meta), path(left),path(right)
    val(bwa_index)

    output:
	tuple val(meta), path(bam), emit: bam
	val(sample), emit: sample_name
    val(meta), emit: meta_data
	path("versions.yml"), emit: versions
    
	script:
	bam = "${meta.sample_id}_${meta.library_id}_${meta.readgroup_id}_bwa-aligned_fm.bam"
	sample = meta.sample_id

	"""
    bwa mem -H ${params.dict} -M -R "@RG\\tID:${meta.readgroup_id}\\tPL:ILLUMINA\\tPU:${meta.platform_unit}\\tSM:${meta.sample_id}\\tLB:${meta.library_id}\\tDS:${bwa_index}\\tCN:${meta.center}" \
        -t ${task.cpus} ${bwa_index} $left $right \
        | samtools fixmate -@ ${task.cpus} -m - - \
		| samtools sort -@ ${task.cpus} -m 4G -O bam -o $bam - 
	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
		
	"""	
}
