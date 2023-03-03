process CNVKIT_EXPORT {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVkit", mode: 'copy'
    
    container 'quay.io/biocontainers/cnvkit:0.9.10--pyhdfd78af_0'

    input:
    tuple val(meta),path(cns)

    output:
    tuple val(meta),path(bed), emit: bed
    tuple val(meta),path(vcf), emit: vcf
    
    script:
    bed = cns.getBaseName() + ".bed"
    vcf = cns.getBaseName() + ".vcf"

    """
        cnvkit.py export bed $cns -y -o $bed
        cnvkit.py export vcf $cns -y -o $vcf
    """

}