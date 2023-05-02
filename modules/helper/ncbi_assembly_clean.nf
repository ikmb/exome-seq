process NCBI_ASSEMBLY_CLEAN {

    tag "${fa}"

    publishDir "${params.outdir}/${assembly}", mode: 'copy'

    container 'quay.io/biocontainers/ruby-dna-tools:1.0--hdfd78af_3'

    input:
    tuple val(meta),path(fa),path(report)

    output:
    tuple val(meta),path(fasta), emit: fasta
    path(look), emit: lookup

    script:
    assembly = meta.assembly
    fasta = assembly + ".fasta"
    look = assembly + "_lookup.txt"
    """
        ncbi_assembly_clean.rb --fasta $fa --report $report --base $assembly

    """

}
