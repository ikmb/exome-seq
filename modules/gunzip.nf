process GUNZIP {

    tag "${gz}"

    input:
    path(gz)

    output:
    path(ungz), emit: decompressed

    script:
    ungz = gz.getBaseName()

    """
        gunzip -c $gz > $ungz
    """

}