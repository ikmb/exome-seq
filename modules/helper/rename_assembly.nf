process RENAME_ASSEMBLY {

    tag "${meta.id}"

    publishDir "${params.outdir}/${meta.assembly}/", mode: 'copy'

    input:
    tuple val(meta),path(fasta)

    output:
    tuple val(meta),path(fasta_n), emit: fasta

    script:
    fasta_n = meta.assembly + ".fasta"

    """
    sed 's/ .*//' $fasta > $fasta_n
    """

}
