process INTERVAL_TO_BED {

	tag "${intervals}"

        executor 'local'

	container 'quay.io/biocontainers/picard:2.20.4--1'

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

