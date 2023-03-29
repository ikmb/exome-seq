process GATK_MAKESITESONLYVCF {

	label 'gatk'

	input:
	tuple val(meta),path(vcf_file),path(tbi)

	output:
	tuple val(meta),path(vcf_sites_only),path(vcf_sites_only_tbi), emit: vcf
	path("versions.yml"), emit: versions

	script:
	vcf_sites_only = vcf_file.getSimpleName() + "-sites_only.vcf.gz"
	vcf_sites_only_tbi = vcf_sites_only + ".tbi"

    """
    gatk MakeSitesOnlyVcf -I $vcf_file -O $vcf_sites_only
    gatk IndexFeatureFile -I $vcf_sites_only

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
