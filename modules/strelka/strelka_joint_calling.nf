process STRELKA_JOINT_CALLING {

    tag "ALL"

    label 'long_parallel'

    publishDir "${params.outdir}/MergedCallset/STRELKA_JOINT_CALLING", mode: 'copy'

    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'

    input:
    tuple path(bams),path(bais)
    tuple path(bed),path(bed_tbi)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple path(vcf_merged),path(vcf_merged_tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    run_dir = "run_dir"
    vcf_merged = "strelka-joint_calling.vcf.gz"
    vcf_merged_tbi = vcf_merged + ".tbi"

    """
    configureStrelkaGermlineWorkflow.py \
        --bam ${bams.join(' --bam ')} \
        --referenceFasta ${fasta} \
        --runDir $run_dir \
        --callRegions $bed \
        --exome

    $run_dir/runWorkflow.py -m local -j ${task.cpus}

    cp $run_dir/results/variants/variants.vcf.gz $vcf_merged
    cp $run_dir/results/variants/variants.vcf.gz.tbi $vcf_merged_tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS

     """

}


