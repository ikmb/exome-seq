process VEP_INSTALL_PLUGINS {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(archive)

    output:
    path("vep")

    script:
    base = archive.getSimpleName()

    """
    tar -xvf $archive
    mv VEP_plugins-release-* vep/plugins
    """
}
