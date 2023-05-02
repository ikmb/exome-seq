process GUNZIP {

    tag "${gz}"

    input:
    tuple val(meta),path(gz)

    output:
    tuple val(meta),path(ungz), emit: decompressed

    script:
    ungz = gz.getBaseName()

    """
        gunzip -c $gz > $ungz
    """

}
