process VEP_INSTALL_CACHE {

    tag "${base}"

    stageOutMode 'rsync' 

    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    path(archive)

    output:
    path("homo_sapiens/*")

    script:
    base = archive.getSimpleName()

    """
    tar -xvf $archive
    """
}
