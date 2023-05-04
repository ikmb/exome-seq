process INTERVAL_TO_BED {

    tag "${intervals}"

    label 'short_serial'

    executor 'local'

    container 'quay.io/biocontainers/picard:2.20.4--1'

    input:
    path(intervals)

    output:
    path(bed_file), emit: bed
    path("versions.yml"), emit: versions

    script:
    bed_file = intervals.getBaseName() + ".bed"

    """
    picard IntervalListTools I=$intervals O=targets.padded.interval_list
    picard IntervalListToBed I=targets.padded.interval_list O=$bed_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: 3.0.0
    END_VERSIONS
    """
}

