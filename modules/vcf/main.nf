process merge_vcf {

	input:
	tuple val(caller),path(vcfs)

	output:
	tuple path(merged_vcf),path(merged_tbi), emit: vcf

	script:
	merged_vcf = caller + ".flat_merged." + params.run_name + ".vcf.gz"
	merged_tbi = merged_vcf + ".tbi"

	"""
		bcftools merge --threads ${task.cpus} -o $merged_vcf -O z *.vcf.gz
		bcftools index -t $merged_vcf
	"""

}

process vcf_get_sample {

        label 'gatk'

        input:
        tuple val(meta),path(vcf),path(vcf_index)
        val(sample_name)

        output:
        tuple path(vcf_sample),path(vcf_sample_index), emit: vcf

        script:
        vcf_sample = sample_name + ".vcf.gz"
        vcf_sample_index = vcf_sample + ".tbi"

        """
                gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $sample_name -O variants.vcf.gz -OVI
                gatk LeftAlignAndTrimVariants -R $FASTA -V variants.vcf.gz -O $vcf_sample
                rm variants.vcf.gz

        """

}

process vcf_add_header {

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Variants", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_r),path(tbi_r), emit: vcf

        script:

        vcf_r = vcf.getBaseName() + ".final.vcf.gz"
        tbi_r = vcf_r + ".tbi"

        """
                echo "##reference=${params.assembly}" > header.txt
                bcftools annotate -h header.txt -O z -o $vcf_r $vcf
                tabix $vcf_r
        """

}

process vcf_stats {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Variants", mode: 'copy'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        path(vcf_stats), emit: stats

        script:
        vcf_stats = vcf.getBaseName() + ".stats"

        """
                bcftools stats $vcf > $vcf_stats
        """

}

process vcf_filter_pass {

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta.patient_id),val(meta.sample_id),path(vcf_pass), path(vcf_pass_index), emit: vcf

	script:
	vcf_pass = vcf.getSimpleName() + ".pass.vcf.gz"
	vcf_pass_index = vcf_pass + ".tbi"

	"""
		bcftools view -f PASS -O z -o $vcf_pass $vcf
		tabix $vcf_pass
	"""

}

process vcf_add_dbsnp {

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_annotated), path(vcf_annotated_index), emit: vcf

        script:
        vcf_annotated = vcf.getBaseName() + ".rsids.vcf.gz"
        vcf_annotated_index = vcf_annotated + ".tbi"

        """
                bcftools annotate -c ID -a $params.dbsnp -O z -o $vcf_annotated $vcf
                tabix $vcf_annotated
        """
}

process vcf_index {

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

process vcf_compress_and_index {

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
