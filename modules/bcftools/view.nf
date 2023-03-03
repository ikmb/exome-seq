process BCFTOOLS_VIEW {

	container 'quay.io/biocontainers/bcftools:1.14--hde04aa1_1'

        tag "${meta.patient_id}|${meta.sample_id}"

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_pass), path(vcf_pass_index), emit: vcf

        script:
        vcf_pass = vcf.getBaseName() + "_pass.vcf.gz"
        vcf_pass_index = vcf_pass + ".tbi"

        """
                bcftools view -f PASS -O z -o $vcf_pass $vcf
                tabix $vcf_pass
        """

}

