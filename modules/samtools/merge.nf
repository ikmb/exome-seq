process SAMTOOLS_MERGE {

        tag "${meta.patient_id}|${meta.sample_id}"

        input:
        tuple val(meta), path(aligned_bam_list)

        output:
        tuple val(meta),path(merged_bam), emit: bam
	val(meta), emit: meta_data

        script:
        merged_bam = meta.patient_id + "_" + meta.sample_id + ".merged.bam"
        merged_bam_index = merged_bam + ".bai"

        """
                        samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}
        """
}

