process VEP_INSTALL_PLUGINS {

    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    path(archive)

    output:
    path("plugins")

    script:
    base = archive.getSimpleName()

    """
    unzip $archive
    mv VEP_plugins-release-* plugins
    """
}
