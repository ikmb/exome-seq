process GATK_GENOMICSDBIMPORT {

	tag "ALL"
	
	label 'gatk'

	input:
	path(gvcfs)
	path(tbis)

	output:
	path(db), emit: db
	path("versions.yml"), emit: versions

	script:
	db = "genomicsdb-" + params.run_name

	"""
	gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
		--merge-input-intervals \
		-V ${gvcfs.join(' -V ')} \
		--genomicsdb-workspace-path $db \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS

	"""	
}
