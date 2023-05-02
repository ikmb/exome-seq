process BGZIP_INDEX {

    container 'quay.io/biocontainers/htslib:1.16--h6bc39ce_0'
    executor 'local'

    input:
    path(bed)

    output:
    tuple path(bed_gz),path(bed_gz_tbi), emit: bedgz
    path("versions.yml"), emit: versions

    script:
    bed_gz = bed.getName() + ".gz"
    bed_gz_tbi = bed_gz + ".tbi"

    """
    bgzip $bed
    tabix -p bed $bed_gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htslib: \$( bgzip --version | head -n1 | sed -e "s/bgzip (htslib) //g" )
    END_VERSIONS
    """
}


