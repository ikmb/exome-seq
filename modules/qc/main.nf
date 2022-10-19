// a simple sex check looking at coverage of the SRY gene
process SEX_CHECK {

        input:
        path(bams)

        output:
        path(sex_check_yaml), emit: yaml

        script:
        sex_check_yaml = "sex_check_mqc.yaml"

        """
                parse_sry_coverage.pl --fasta ${params.fasta} --region ${params.sry_region} > $sex_check_yaml
        """
}

process QC_SNPS {

	input:
	tuple val(meta),path(vcf),path(tbi)
	path(bed)

	output:
	tuple val(meta),path(json)

	script:
	
	json = vcf.getBaseName() + ".json"

	"""
		bcftools view -T $bed $vcf > snps.vcf
	"""
}
