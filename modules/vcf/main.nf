process MERGE_VCF {

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

process VCF_GATK_SORT {

	label 'picard'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_sorted),path(tbi_sorted), emit: vcf

	script:
	vcf_sorted = vcf.getSimpleName() + ".sorted.vcf.gz"
	tbi_sorted = vcf_sorted + ".tbi"

	"""
		picard SortVcf I=$vcf O=$vcf_sorted CREATE_INDEX=true 
		
	"""

}
	
process VCF_GET_SAMPLE {

        label 'gatk'

        input:
        tuple val(m_f),path(vcf),path(vcf_index)
        val(meta)

        output:
        tuple val(meta),path(vcf_sample),path(vcf_sample_index), emit: vcf

        script:
        def prefix = meta.patient_id + "_" + meta.sample_id
        vcf_sample = prefix + "-" + m_f.variantcaller + ".split.vcf.gz"
        vcf_sample_index = vcf_sample + ".tbi"

        """
                gatk SelectVariants --remove-unused-alternates --exclude-non-variants -V $vcf -sn $prefix -O $vcf_sample -OVI

        """

}

process VCF_ADD_HEADER {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

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

process VCF_STATS {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Stats", mode: 'copy'

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

process VCF_FILTER_PASS {

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_pass), path(vcf_pass_index), emit: vcf

	script:
	vcf_pass = vcf.getSimpleName() + ".pass.vcf.gz"
	vcf_pass_index = vcf_pass + ".tbi"

	"""
		bcftools view -f PASS -O z -o $vcf_pass $vcf
		tabix $vcf_pass
	"""

}

process VCF_ADD_DBSNP {

	//publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Variants", mode: 'copy'

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

process VCF_INDEX {

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

process VCF_COMPRESS_AND_INDEX {

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
