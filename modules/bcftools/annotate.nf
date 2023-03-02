process BCFTOOLS_ANNOTATE {

	container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

        tag "${meta.patient_id}|${meta.sample_id}"

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_r),path(tbi_r), emit: vcf

        script:

        vcf_r = vcf.getBaseName() + "-final.vcf.gz"
        tbi_r = vcf_r + ".tbi"

        """
                echo "##reference=${params.assembly}" > header.txt
                bcftools annotate -h header.txt -O z -o $vcf_r $vcf
                tabix $vcf_r
        """

}


