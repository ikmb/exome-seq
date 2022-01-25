process expansion_hunter {

	label 'expansion_hunter'

	publishDir "${params.outdir}/${indivID}/${sampleID}/ExpansionHunter", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID), path(bam),path(bai)
	path(catalog)

	output:
	tuple val(indivID), val(sampleID),path(expansion_report)
	path(expansion_vcf)

	script:
	expansion_report = indivID + "_" + sampleID + ".expansion_report.json"
	expansion_vcf = indivID + "_" + sampleID + ".expansion_report.vcf"
	prefix = indivID + "_" + sampleID + ".expansion_report"

	"""
		ExpansionHunter --reads $bam --reference $params.fasta --variant-catalog $catalog --output-prefix $prefix
	"""

}

process  expansions2xls {

	publishDir "${params.outdir}/${indivID}/${sampleID}/ExpansionHunter", mode: 'copy'

	input:
	tuple val(indivID), val(sampleID),path(report)

	output:
	path(expansion_xls)

	script:

	expansion_xls = report.getBaseName() + ".xlsx"

	"""
		expansions2xls.pl --infile $report --outfile $expansion_xls
	"""

}

