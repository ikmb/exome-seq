process BEDTOOLS_SLOP {

    tag "${b}"

    container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

    input:
    path(b)
    val(padding)
    tuple path(fasta),path(fai),path(dict)

    output:
    path(bed_padded), emit: bed
    path("versions.yml"), emit: versions

    script:
    bed_padded = b.getBaseName() + ".padded.bed"

    """
    bedtools slop -i $b -g $fai -b $padding > $bed_padded

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$( bedtools -version | head -n1 | sed -e "s/bedtools v//g" )
    END_VERSIONS
    """
}


