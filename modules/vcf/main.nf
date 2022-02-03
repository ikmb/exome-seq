process merge_vcf {

	input:
	path(vcfs)

	output:
	tuple val("Deepvariant"),val("Merged"),val("Bcftools"),path(merged_vcf),path(merged_tbi)

	script:
	merged_vcf = "deepvariant.flat_merged." + params.run_name + ".vcf.gz"
	merged_tbi = merged_vcf + ".tbi"

	"""
		bcftools merge --threads ${task.cpus} -o $merged_vcf -O z *.vcf.gz
		bcftools index -t $merged_vcf
	"""

}

process vcf_get_sample {

        label 'gatk'

        input:
        tuple path(vcf),path(vcf_index)
        val(sample_name)

        output:
        set file(vcf_sample),file(vcf_sample_index)

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

        publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'

        input:
        tuple val(cname),val(indivID),val(sampleID),path(vcf),path(tbi)

        output:
        tuple val(cname),val(indivID),val(sampleID),path(vcf_r),path(tbi_r)

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

	publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'

        input:
        tuple val(cname),val(indivID),val(sampleID),path(vcf),path(tbi)

        output:
        path(vcf_stats)

        script:
        vcf_stats = vcf.getBaseName() + ".stats"

        """
                bcftools stats $vcf > $vcf_stats
        """

}

process vcf_filter_pass {

	input:
	tuple val(cname),val(indivID),val(sampleID),path(vcf),path(tbi)

	output:
	tuple val(cname),val(indivID),val(sampleID),path(vcf_pass), path(vcf_pass_index)

	script:
	vcf_pass = vcf.getSimpleName() + ".pass.vcf.gz"
	vcf_pass_index = vcf_pass + ".tbi"

	"""
		bcftools view -f PASS -O z -o $vcf_pass $vcf
		tabix $vcf_pass
	"""

}

process vcf_add_dbsnp {

        //publishDir "${params.outdir}/${indivID}/${sampleID}/Variants", mode: 'copy'

        input:
        tuple val(cname),val(indivID),val(sampleID),path(vcf),path(tbi)

        output:
        tuple val(cname),val(indivID),val(sampleID),path(vcf_annotated), path(vcf_annotated_index)

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
	tuple val(cname),val(indivID),val(sampleID),path(vcf)

	output:
	tuple val(cname),val(indivID),val(sampleID),path(vcf),path(tbi)

	script:

	tbi = vcf + ".tbi"

	"""
		tabix $vcf
	"""

}

process vcf_compress_and_index {

	input:
	tuple val(indivID),val(sampleID),path(vcf)

	output:
	tuple val(indivID),val(sampleID),path(vcf_gz),path(vcf_gz_tbi)

	script:
	vcf_gz = vcf + ".gz"
	vcf_gz_tbi = vcf_gz + ".tbi"

	"""
		bgzip -c $vcf > $vcf_gz
		tabix $vcf_gz
	"""
}
