process GLNEXUS {

	container 'quay.io/mlin/glnexus:v1.3.1'

        publishDir "${params.outdir}/MergedCallset/GLNEXUS_DEEPVARIANT", mode: 'copy'

        input:
        path(gvcfs)
        path(bed)

        output:
        tuple path(merged_vcf),path(merged_tbi), emit: vcf

        script:
        merged_vcf = "deepvariant.joint_merged." + params.run_name + ".vcf.gz"
        merged_tbi = merged_vcf + ".tbi"

        """
                /usr/local/bin/glnexus_cli \
                --config ${params.glnexus_config} \
                --bed $bed \
                $gvcfs | bcftools view - | bgzip -c > $merged_vcf
                tabix $merged_vcf

         """
}

