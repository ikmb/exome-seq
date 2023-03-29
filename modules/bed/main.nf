process BED_TO_BEDGZ {

    tag "${bed}"

    executor 'local'

    input:
    path(bed)

    output:
    tuple path(bed_gz),path(bed_gz_tbi), emit: bedgz
    path("versions.yml"), emit: versions

    script:
    bed_gz = bed.getBaseName() + ".bed.gz"
    bed_gz_tbi = bed_gz + ".tbi"

    """
    bgzip $bed
    tabix -p bed $bed_gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$( bgzip --version | head -n1 | sed -e "s/bgzip \(htslib\) //g" )
    END_VERSIONS
    """
}


