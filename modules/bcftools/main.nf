process CONCAT {

	publishDir "${params.outdir}/ConcatenatedVariants", mode: 'copy'

	label 'bcftools'

	input:
	path(vcfs)

	output:
	tuple val(meta), path(merged_vcf),path(merged_vcf_tbi), emit: vcf

	script:
	meta = [:]
	meta.variant_caller = "ConcatenatedCallsets"
	
	merged_vcf = "merged_callset." + params.run_name + ".vcf.gz"
	merged_vcf_tbi = merged_vcf + ".tbi"

	"""
		bcftools concat -a -D -o $merged_vcf -O z *.vcf.gz
		bcftools index -t $merged_vcf
	"""

}

process VCF_SORT {
	
	label 'bcftools'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_sorted),path(tbi_sorted), emit: vcf

	script:
	vcf_sorted = vcf.getSimpleName() + ".sorted.vcf.gz"
	tbi_sorted = vcf_sorted + ".tbi"

	"""
		bcftools sort -o $vcf_sorted $vcf 
		bcftools index -t $vcf_sorted
	"""

}

process VCF_INDEL_NORMALIZE {

        label 'bcftools'

        input:
        tuple val(meta),path(vcf),path(tbi)

        output:
        tuple val(meta),path(vcf_norm),path(tbi_norm), emit: vcf

        script:
        vcf_norm = vcf.getSimpleName() + ".normalized.vcf.gz"
        tbi_norm = vcf_norm + ".tbi"

        """
                bcftools norm -O z -f $params.fasta -o $vcf_norm $vcf
                bcftools index -t $vcf_norm
        """



}

process CSQ {

	label 'bcftools'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/${meta.variantcaller}", mode: 'copy'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple val(meta),path(vcf_fixed),path(tbi_fixed), emit: vcf

	script:
	vcf_fixed = vcf.getSimpleName() + "_csq.vcf.gz"
	tbi_fixed = vcf_fixed + ".tbi"

	"""
		bcftools csq -f $params.fasta -g $params.gtf --phase a $vcf -o $vcf_fixed
		bcftools index -t $vcf_fixed
	"""

}

