process GATK_SELECTVARIANTS {

        tag "${meta.patient_id}|${meta.sample_id}"

        label 'gatk'

        input:
        tuple val(m_f),path(vcf),path(vcf_index)
        val(meta)

        output:
        tuple val(meta),path(vcf_sample),path(vcf_sample_index), emit: vcf

        script:
        def prefix = meta.patient_id + "_" + meta.sample_id
        vcf_sample = prefix + "-" + m_f.variantcaller + "-split.vcf.gz"
        vcf_sample_index = vcf_sample + ".tbi"

        """
                gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $prefix -O $vcf_sample -OVI

        """

}

