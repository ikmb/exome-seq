process HTSLIB_BGZIP {

	container 'quay.io/biocontainers/htslib:1.16--h6bc39ce_0'

        tag "${meta.patient_id}|${meta.sample_id}"

        input:
        tuple val(meta),path(vcf)

        output:
        tuple val(meta),path(vcf_gz),path(vcf_gz_tbi), emit: vcf

        script:
        vcf_gz = vcf + ".gz"
        vcf_gz_tbi = vcf_gz + ".tbi"

        """
                bgzip -c $vcf > $vcf_gz
                tabix $vcf_gz
        """
}

