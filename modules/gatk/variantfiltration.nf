process GATK_VARIANTFILTRATION {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'

	label 'gatk'

	input:
	tuple val(meta),path(vcf),path(tbi)

	output:
	tuple path(vcf_filtered),path(vcf_filtered_tbi), emit: vcf
	path("versions.yml"), emit: versions

	script:
	vcf_filtered = vcf.getBaseName() + "-filtered.vcf.gz"
	vcf_filtered_tbi = vcf_filtered +".tbi"

    """
    gatk --java-options "-Xmx3g -Xms3g" VariantFiltration \
        -V $vcf \
        --filter-expression "ExcessHet > 54.69" \
        --filter-name ExcessHet \
        -O $vcf_filtered -OVI true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
