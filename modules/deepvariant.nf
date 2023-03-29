process DEEPVARIANT {

    tag "${meta.patient_id}|${meta.sample_id}"

    container 'google/deepvariant:1.5.0'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/DEEPVARIANT", mode: 'copy'

    input:
    tuple val(meta), path(bam),path(bai)
    path(bed)
    tuple path(fai),path(fastagz),path(gzfai),path(gzi)

    output:
    path(dv_gvcf), emit: gvcf
    tuple val(meta),path(dv_vcf), emit: vcf
    val(sample_name), emit: sample_name
    path("versions.yml"), emit: versions

    script:
    dv_gvcf = meta.patient_id + "_" + meta.sample_id + "-deepvariant.g.vcf.gz"
    dv_vcf = meta.patient_id + "_" + meta.sample_id + "-deepvariant.vcf.gz"
    sample_name = "${meta.patient_id}_${meta.sample_id}"
    """
    /opt/deepvariant/bin/run_deepvariant \
        --model_type=WES \
        --ref=$fastagz \
        --reads $bam \
        --output_vcf=$dv_vcf \
        --output_gvcf=$dv_gvcf \
        --regions=$bed \
        --num_shards=${task.cpus} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(echo \$(/opt/deepvariant/bin/run_deepvariant --version) | sed 's/^.*version //; s/ .*\$//' )
    END_VERSIONS
    """
}

