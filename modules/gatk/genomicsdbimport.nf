process GATK_GENOMICSDBIMPORT {

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
			-V ${gvcfs.join(' -V ')} \
			--genomicsdb-workspace-path $db \
			-L chr22

	"""	
}