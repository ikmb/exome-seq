process BWA2_MEM {

	tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0'

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
	bam = "${meta.sample_id}_${meta.library_id}_${meta.readgroup_id}_bwa2-aligned.fm.bam"
	sample = "${meta.patient_id}_${meta.sample_id}"

	"""
    bwa-mem2 mem -K 1000000 -H ${params.dict} -M -R "@RG\\tID:${meta.readgroup_id}\\tPL:ILLUMINA\\tPU:${meta.platform_unit}\\tSM:${meta.patient_id}_${meta.sample_id}\\tLB:${meta.library_id}\\tDS:${bwa_index}\\tCN:${meta.center}" \
        -t ${task.cpus} ${bwa_index} $left $right \
        | samtools fixmate -@ ${task.cpus} -m - - \
        | samtools sort -@ ${task.cpus} -m 2G -O bam -o $bam - 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
	"""	
}
