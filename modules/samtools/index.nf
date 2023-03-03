process SAMTOOLS_INDEX {

	container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

        tag "${meta.patient_id}|${meta.sample_id}"

        input:
        tuple val(meta), path(bam)

        output:
        tuple val(meta), path(bam),path(bam_index), emit: bam

        script:
        bam_index = bam.getName() + ".bai"

        """
                samtools index $bam
        """

}

