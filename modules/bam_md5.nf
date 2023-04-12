process BAM_MD5 {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/", mode: 'copy'

    input:
    tuple val(meta),path(bam),path(bai)

    output:
    tuple val(meta),path(bam),path(bai), emit: bam
    path(md5)

    script:
    md5 = bam + ".md5"

    """
        md5sum $bam > $md5
    """
}