process GATK_FILTERVARIANTTRANCHES {

	tag "${meta.patient_id}|${meta.sample_id}"

        container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'medium_serial'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'
	
	input:
	tuple val(meta),path(vcf),path(tbi)
	path(snps)
	path(snps_tbi)
	path(indels)
	path(indels_tbi)

	output:
	tuple val(meta),path(vcf_filtered),path(vcf_filtered_tbi), emit: vcf
	path("versions.yml"), emit: versions

	script:
	
	vcf_filtered = vcf.getSimpleName() + "-filtered_tranches.vcf.gz"
	vcf_filtered_tbi = vcf_filtered + ".tbi"

    """
    gatk FilterVariantTranches \
        -V $vcf \
        --resource ${snps.join(' --resource ')} \
        --resource ${indels.join(' --resource ')} \
        --info-key CNN_1D \
        --invalidate-previous-filters \
        --snp-tranche 99.9 --indel-tranche 99.9 \
        -O $vcf_filtered

    gatk IndexFeatureFile -I $vcf_filtered

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
