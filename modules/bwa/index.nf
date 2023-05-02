process BWA_MEM_INDEX {

    tag "${fa}"

    label 'high_mem'

    publishDir "${params.outdir}/${assembly}/bwa", mode: 'copy'

    input:
    tuple val(meta), path(fa)

    output:
    tuple path(fa),path("*.ann"),path("*.bwt"),path("*.pac"),path("*.sa"), emit: bwa_index
    path("versions.yml"), emit: versions
    
    script:
    assembly = meta.assembly
    """
    bwa index $fa
	
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
		
	"""	
}
