process GATK_GENOMICSDBIMPORT {

	tag "ALL"
	
	label 'gatk'

	input:
	path(gvcfs)
	path(tbis)

	output:
	path(db), emit: db

	script:
	db = "genomicsdb-" + params.run_name

	"""
		gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
			--merge-input-intervals \
			-V ${gvcfs.join(' -V ')} \
			--genomicsdb-workspace-path $db \

	"""	
}
