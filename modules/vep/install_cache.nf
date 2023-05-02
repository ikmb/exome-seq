process VEP_INSTALL_CACHE {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(archive)

    output:
    path("vep")

    script:
    base = archive.getSimpleName()

    """
    tar -xvf $archive
    mv homo_sapiens vep/homo_sapiens
    """
}
