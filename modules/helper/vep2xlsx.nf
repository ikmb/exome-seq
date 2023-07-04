process VEP2XLSX {

    container 'ikmb/exome-seq:devel'

    publishDir "${params.outdir}/VEP"
    
    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path(xlsx), emit: xlsx

    script:
    xlsx = vcf.getBaseName() + ".xlsx"

    """
    ruby vep2xlsx.rb -i $vcf -o $xlsx
    """

}