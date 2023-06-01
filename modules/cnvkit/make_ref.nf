process CNVKIT_MAKE_REF {

    label 'medium_parallel'

    publishDir "${params.outdir}/CNVkit", mode: 'copy'

    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    path(bams)
    path(bed)
    tuple path(fasta),path(fai),path(dict)
    val(nsamples)
   
    output:
    path(cnn_f), emit: cnn

    script:
   
    cnn_f = params.run_name + "-cnvkit.cnn"

    if (nsamples >= 5) {
        """
        cnvkit.py access $fasta -x ${baseDir}/assets/cnvkit/GRCh38/hg38.encode_exclusion_list.bed -o access.bed
        cnvkit.py antitarget $bed -g access.bed -o antitargets.bed
        cnvkit.py batch -n *.bam -o $cnn_f --targets $bed -a antitargets.bed -f $fasta -d results -p ${task.cpus}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            cnvkit: \$( cnvkit.py version )
        END_VERSIONS

        """
    } else {

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
}
