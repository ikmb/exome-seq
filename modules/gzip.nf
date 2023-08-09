process GZIP {

    tag "${f}"

    input:
    tuple val(meta),path(f)

    output:
    tuple val(meta),path(gz), emit: compressed

    script:
    gz = f.getName() + ".gz"

    """
        gzip -c $f > $gz
    """

}
