process BCFTOOLS_CONCAT {

        label 'bcftools'

        publishDir "${params.outdir}/ConcatenatedVariants", mode: 'copy'

        label 'bcftools'

        input:
        path(vcfs)

        output:
        tuple val(meta), path(merged_vcf),path(merged_vcf_tbi), emit: vcf

        script:
        meta = [:]
        meta.variant_caller = "ConcatenatedCallsets"

        merged_vcf = "merged_callset." + params.run_name + ".vcf.gz"
        merged_vcf_tbi = merged_vcf + ".tbi"

        """
                bcftools concat -a -D -o $merged_vcf -O z *.vcf.gz
                bcftools index -t $merged_vcf
        """

}

