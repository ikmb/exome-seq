process CNVKIT_MAKE_FLAT_REF {

    label 'medium_parallel'

    publishDir "${params.outdir}/CNVkit", mode: 'copy'

    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    path(bed)
    tuple path(fasta),path(fai),path(dict)
   
    output:
    path(cnn_f), emit: cnn

    script:
   
    cnn_f = params.run_name + "-flatref-cnvkit.cnn"

    """
    cnvkit.py access $fasta -x ${baseDir}/assets/cnvkit/GRCh38/hg38.encode_exclusion_list.bed -o access.bed
    cnvkit.py antitarget $bed -g access.bed -o antitargets.bed
    cnvkit.py reference -o $cnn_f -f $fasta -t $bed -a antitargets.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$( cnvkit.py version )
    END_VERSIONS
    """
}
