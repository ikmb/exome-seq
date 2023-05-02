process DRAGMAP_INDEX {

    tag "${fa}"

    publishDir "${params.outdir}/${assembly}", mode: 'copy'

    label 'high_mem'

    container 'quay.io/biocontainers/mulled-v2-580d344d9d4a496cd403932da8765f9e0187774d:5ebebbc128cd624282eaa37d2c7fe01505a91a69-0'

    input:
    tuple val(meta),path(fa)

    output:
    path("dragmap"), emit: dragen_index
    path("versions.yml"), emit: versions

    script:	
    assembly = meta.assembly
    """
    mkdir -p dragmap

    dragen-os \
        --build-hash-table true \
        --ht-reference $fa \
        --output-directory dragmap \
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragmap: \$(echo \$(dragen-os --version 2>&1))
    END_VERSIONS

	"""
}
