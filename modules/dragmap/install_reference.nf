process DRAGMAP_INSTALL_REFERENCE {

    publishDir "${params.outdir}/${meta.assembly}", mode: 'copy'

    input:
    tuple val(meta),path(archive)

    output:
    tuple val(meta),path("dragmap"), emit: ref_dir

    script:

    """
    mkdir -p dragmap
    cp $archive dragmap/
    cd dragmap
    bash $archive
    """
}