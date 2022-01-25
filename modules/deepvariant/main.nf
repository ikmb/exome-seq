process deepvariant {

        label 'deepvariant'

        publishDir "${params.outdir}/${indivID}/${sampleID}/DeepVariant", mode: 'copy'

        input:
        tuple val(indivID), val(sampleID), path(bam),path(bai)
        path(bed)
	tuple path(fai),path(fastagz),path(gzfai),path(gzi)

        output:
        path(gvcf)
        tuple val("Deepvariant"),val(indivID),val(sampleID),path(vcf)
	path(vcf)
        val(sample_name)

        script:
        gvcf = bam.getBaseName() + ".g.vcf.gz"
        vcf = bam.getBaseName() + ".vcf.gz"
        sample_name = indivID + "_" + sampleID

        """
                /opt/deepvariant/bin/run_deepvariant \
                --model_type=WES \
                --ref=$fastagz \
                --reads $bam \
                --output_vcf=$vcf \
                --output_gvcf=$gvcf \
                --regions=$bed \
                --num_shards=${task.cpus} \
        """
}


process merge_gvcfs {

	label 'glnexus'

	input:
	path(gvcfs)
	bed(bed)

	output:
	tuple val("Deepvariant"),val("JointCalling"),val("GLNexus"),path(merged_vcf)

	script:
	merged_vcf = "deepvariant.joint_merged." + run_name + ".vcf.gz"

	"""
		/usr/local/bin/glnexus_cli \
		--config ${params.glnexus_config} \
		--bed $bed \
		$gvcfs | bcftools view - | bgzip -c > $merged_vcf

	 """
}
