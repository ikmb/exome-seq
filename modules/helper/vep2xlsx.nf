process VEP2XLSX {

    container 'ikmb/exome-seq:devel'

    publishDir "${params.outdir}/VEP", mode: 'copy'
    
    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(sheet), emit: xlsx

    script:
    sheet = vcf.getBaseName() + ".xlsx"

    """
    vep2xls_fast.rb -i $vcf -o $sheet
    """

}
