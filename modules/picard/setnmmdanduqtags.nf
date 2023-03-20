process PICARD_SET_BAM_TAGS {

        tag "${meta.patient_id}|${meta.sample_id}"

	container 'quay.io/biocontainers/picard:3.0.0--hdfd78af_1'
	
        input:
        tuple val(meta),path(bam),path(bai)
	tuple path(fasta),path(fai),path(dict)

        output:
        tuple val(meta),path(bam_fixed), emit: bam

        script:
        bam_fixed = bam.getBaseName() + "-fixed.bam"
        bai_fixed = bam_fixed + ".bai"

        """
        picard SetNmMdAndUqTags -Xmx${task.memory.toGiga()}G \
                R=$fasta \
                I=$bam \
                O=$bam_fixed
        """
}

