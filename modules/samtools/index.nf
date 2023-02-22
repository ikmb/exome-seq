process SAMTOOLS_INDEX {

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

