process PICARD_SORTVCF {

        label 'picard'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_sorted),path(tbi_sorted), emit: vcf

        script:
        vcf_sorted = vcf.getSimpleName() + ".sorted.vcf.gz"
        tbi_sorted = vcf_sorted + ".tbi"

        """
                picard SortVcf I=$vcf O=$vcf_sorted CREATE_INDEX=true

        """

}

