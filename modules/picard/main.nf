process interval_to_bed {

        executor 'local'

        label 'picard'

        input:
        path(intervals)

        output:
        path(bed)

        script:
        bed = intervals.getBaseName() + ".bed"

        """
                picard IntervalListTools I=$intervals O=targets.padded.interval_list PADDING=$params.interval_padding
                picard IntervalListToBed I=targets.padded.interval_list O=$bed
        """
}

process multi_metrics {

        label 'picard'

        publishDir "${params.outdir}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(bam), path(bai)
        path(baits)

        output:
        file("${prefix}*") 

        script:
        prefix = indivID + "_" + sampleID + "."

        """
                picard -Xmx5g CollectMultipleMetrics \
                PROGRAM=MeanQualityByCycle \
                PROGRAM=QualityScoreDistribution \
                PROGRAM=CollectAlignmentSummaryMetrics \
                PROGRAM=CollectInsertSizeMetrics\
                PROGRAM=CollectSequencingArtifactMetrics \
                PROGRAM=CollectQualityYieldMetrics \
                PROGRAM=CollectGcBiasMetrics \
                PROGRAM=CollectBaseDistributionByCycle \
                INPUT=${bam} \
                REFERENCE_SEQUENCE=${params.fasta} \
                DB_SNP=${params.dbsnp} \
                INTERVALS=${baits} \
                ASSUME_SORTED=true \
                QUIET=true \
                OUTPUT=${prefix} \
                TMP_DIR=tmp
        """
}

process hybrid_capture_metrics {

        label 'picard'

        publishDir "${params.outdir}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(bam), path(bai)
        path(targets)
        path(baits)

        output:
        path(outfile)
        path(outfile_per_target)

        script:
        outfile = indivID + "_" + sampleID + ".hybrid_selection_metrics.txt"
        outfile_per_target = indivID + "_" + sampleID + ".hybrid_selection_per_target_metrics.txt"

        """
        picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                PER_TARGET_COVERAGE=${outfile_per_target} \
                TARGET_INTERVALS=${targets} \
                BAIT_INTERVALS=${baits} \
                REFERENCE_SEQUENCE=${params.fasta} \
                MINIMUM_MAPPING_QUALITY=$params.min_mapq \
                TMP_DIR=tmp
        """
}

process oxo_metrics {


	label 'picard'

	publishDir "${params.outdir}/${indivID}/${sampleID}/Processing/Picard_Metrics", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID), path(bam), path(bai)
	path(targets)

	output:
	path(outfile)

	script:
	outfile = indivID + "_" + sampleID + ".OxoG_metrics.txt"

	"""

         picard -Xmx${task.memory.toGiga()-1}G CollectOxoGMetrics \
                INPUT=${bam} \
                OUTPUT=${outfile} \
                DB_SNP=${params.dbsnp} \
                INTERVALS=${targets} \
                REFERENCE_SEQUENCE=${params.fasta} \
                TMP_DIR=tmp
        """
}

process panel_coverage {

        publishDir "${params.outdir}/Summary/Panel/PanelCoverage", mode: "copy"

        input:
        tuple val(indivID),val(sampleID),file(bam),file(bai),file(panel)
        path(targets)

        output:
        tuple val(panel_name),path(coverage)
        tuple val(indivID),val(sampleID),path(target_coverage_xls)
        path(target_coverage)

        script:
        panel_name = panel.getSimpleName()
        coverage = indivID + "_" + sampleID + "." +  panel_name  + ".hs_metrics.txt"
        target_coverage = indivID + "_" + sampleID + "." +  panel_name  + ".per_target.hs_metrics.txt"
        target_coverage_xls = indivID + "_" + sampleID + "." + panel_name + ".per_target.hs_metrics_mqc.xlsx"

        // optionally support a kill list of known bad exons
        def options = ""
        if (params.kill) {
                options = "--ban ${params.kill}"
        }

        // do something here - get coverage and build an XLS sheet
        // First we identify which analysed exons are actually part of the exome kit target definition.
        """

                picard -Xmx${task.memory.toGiga()}G IntervalListTools \
                        INPUT=$panel \
                        SECOND_INPUT=$targets \
                        ACTION=SUBTRACT \
                        OUTPUT=overlaps.interval_list

                picard -Xmx${task.memory.toGiga()}G CollectHsMetrics \
                        INPUT=${bam} \
                        OUTPUT=${coverage} \
                        TARGET_INTERVALS=${panel} \
                        BAIT_INTERVALS=${panel} \
                        CLIP_OVERLAPPING_READS=false \
                        REFERENCE_SEQUENCE=${params.fasta} \
                        TMP_DIR=tmp \
                        MINIMUM_MAPPING_QUALITY=$params.min_mapq \
                        PER_TARGET_COVERAGE=$target_coverage

                target_coverage2xls.pl $options --infile $target_coverage --min_cov $params.panel_coverage --skip overlaps.interval_list --outfile $target_coverage_xls

        """
}
