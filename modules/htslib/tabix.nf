process TABIX {

	container 'quay.io/biocontainers/htslib:1.16--h6bc39ce_0'

        tag "${meta.patient_id}|${meta.sample_id}"

        input:
        tuple val(meta),path(vcf)

        output:
        tuple val(meta),path(vcf),path(tbi), emit: vcf

        script:

        tbi = vcf + ".tbi"

        """
                tabix $vcf
        """

}

