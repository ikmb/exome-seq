process DEEPVARIANT {

        label 'deepvariant'

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/DEEPVARIANT", mode: 'copy'

        input:
        tuple val(meta), path(bam),path(bai)
        path(bed)
	tuple path(fai),path(fastagz),path(gzfai),path(gzi)

        output:
        path(dv_gvcf), emit: gvcf
        tuple val(meta),path(dv_vcf), emit: vcf
        val(sample_name), emit: sample_name

        script:
        dv_gvcf = meta.patient_id + "_" + meta.sample_id + "-deepvariant.g.vcf.gz"
        dv_vcf = meta.patient_id + "_" + meta.sample_id + "-deepvariant.vcf.gz"
        sample_name = "${meta.patient_id}_${meta.sample_id}"
        """
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=WES \
                --ref=$fastagz \
                --reads $bam \
                --output_vcf=$dv_vcf \
                --output_gvcf=$dv_gvcf \
                --regions=$bed \
                --num_shards=${task.cpus} \
        """
}


process MERGE_GVCFS {

	publishDir "${params.outdir}/MergedCallset/GLNEXUS_DEEPVARIANT", mode: 'copy'

	label 'glnexus'

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
