process GATK_MERGEVCFS {

    tag "${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MergeVCF", mode: 'copy'    

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    input:
    tuple val(meta),path(vcfs),path(tbis)

    output:
    tuple val(meta),path(merged_vcf),path(merged_vcf_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    merged_vcf = meta.sample_id + "-merged.vcf.gz"
    merged_vcf_tbi = merged_vcf + ".tbi"

    """
    gatk --java-options "-Xmx4g" MergeVcfs \
        --INPUT ${vcfs.sort().join(' --INPUT ')} \
        --OUTPUT $merged_vcf

    gatk IndexFeatureFile -I $merged_vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
