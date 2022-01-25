process merge_multi_lane {

        input:
        tuple val(indivID), val(sampleID), path(aligned_bam_list)

        output:
        tuple val(indivID),val(sampleID),path(merged_bam)

        script:
        merged_bam = indivID + "_" + sampleID + ".merged.bam"
        merged_bam_index = merged_bam + ".bai"

        """
                        samtools merge -@ 4 $merged_bam ${aligned_bam_list.join(' ')}
        """
}

process bam_index {

        publishDir "${params.outdir}/${indivID}/${sampleID}/", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(bam)

        output:
        tuple val(indivID), val(sampleID), path(bam),path(bam_index)
        tuple path(bam),path(bam_index)

        script:
        bam_index = bam.getName() + ".bai"

        """
                samtools index $bam
        """

}

process dedup {

        publishDir "${params.outdir}/${indivID}/${sampleID}/", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(merged_bam),path(merged_bam_index)

        output:
        tuple val(indivID), val(sampleID), path(outfile_bam),path(outfile_bai)
        tuple path(outfile_bam),path(outfile_bai)
        path(outfile_md5)
        path(outfile_metrics)

        script:
        outfile_bam = indivID + "_" + sampleID + ".dedup.bam"
        outfile_bai = indivID + "_" + sampleID + ".dedup.bam.bai"
        outfile_md5 = indivID + "_" + sampleID + ".dedup.bam.md5"
        outfile_metrics = indivID + "_" + sampleID + "_duplicate_metrics.txt"

        """
                samtools markdup -@ ${task.cpus} $merged_bam $outfile_bam
                samtools index $outfile_bam
                samtools stats $outfile_bam > $outfile_metrics
                md5sum $outfile_bam > $outfile_md5
        """

}

