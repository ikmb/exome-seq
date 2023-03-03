process BCFTOOLS_MERGE {

	container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

        publishDir "${params.outdir}/MergedCallset/Bcftools", mode: 'copy'

        input:
        tuple val(meta),path(vcfs),path(tbis)

        output:
        tuple val(meta),path(merged_vcf),path(merged_tbi), emit: vcf

        script:
        merged_vcf = meta.variantcaller + ".flat_merged." + params.run_name + ".vcf.gz"
        merged_tbi = merged_vcf + ".tbi"
        """
                bcftools merge --threads ${task.cpus} -o $merged_vcf -O z *.vcf.gz
                bcftools index -t $merged_vcf
        """

}

