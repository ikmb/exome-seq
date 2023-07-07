process VEP2XLSX {

    container 'ikmb/exome-seq:5.2'

    publishDir "${params.outdir}/VEP", mode: 'copy'
    
    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(sheet), emit: xlsx

    script:
    sheet = vcf.getBaseName() + ".xlsx"

    """
    vep2xls.rb -i $vcf -o $sheet
    """

}
