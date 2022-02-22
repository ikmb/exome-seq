process merge_multi_lane {

        input:
        tuple val(meta), path(aligned_bam_list)

        output:
        tuple val(meta),path(merged_bam), emit: bam

        script:
        merged_bam = meta.patient_id + "_" + meta.sample_id + ".merged.bam"
        merged_bam_index = merged_bam + ".bai"

        """
                        samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}
        """
}

process bam_index {

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

process dedup {

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/", mode: 'copy'

        input:
        tuple val(meta),path(merged_bam),path(merged_bam_index)

        output:
        tuple val(meta), path(outfile_bam),path(outfile_bai), emit: bam
        path(outfile_md5), emit: md5sum
        path(outfile_metrics), emit: report

        script:
        prefix = meta.patient_id + "_" + meta.sample_id + ".dedup"
        outfile_bam = prefix + ".bam"
        outfile_bai = prefix + ".bam.bai"
        outfile_md5 = prefix + ".bam.md5"
        outfile_metrics = prefix + "_duplicate_metrics.txt"

        """
                samtools markdup -@ ${task.cpus} $merged_bam $outfile_bam
                samtools index $outfile_bam
                samtools stats $outfile_bam > $outfile_metrics
                md5sum $outfile_bam > $outfile_md5
        """

}
