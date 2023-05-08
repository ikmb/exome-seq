process VEP_INSTALL_CACHE {

    executor 'local'

    tag "${base}"

    label 'medium_serial'

    stageOutMode 'rsync' 

    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    val(archive)

    output:
    path("homo_sapiens/*")

    script:
    base = file(archive).getName()

    """
    wget $archive
    tar -xvf $base
    """
}
