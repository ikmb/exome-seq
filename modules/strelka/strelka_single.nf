process STRELKA_SINGLE_SAMPLE {

    tag "${meta.patient_id}|${meta.sample_id}"

    label 'long_parallel'

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/STRELKA", mode: 'copy'

    container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'

    input:
    tuple val(meta),path(bam),path(bai),path(indels),path(indels_tbi)
    tuple path(bed),path(bed_tbi)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(vcf), emit: vcf
	path("versions.yml"), emit: versions

	script:
    def options = ""
    if (indels) {
        options = "--indelCandidates $indels"
    }
    run_dir = "run_dir"
    vcf = meta.patient_id + "_" + meta.sample_id + "-strelka.vcf.gz"
    tbi = vcf + ".tbi"

    """
    configureStrelkaGermlineWorkflow.py \
        --bam $bam \
        --referenceFasta ${fasta} \
        --runDir $run_dir \
        --callRegions $bed \
        --exome \
        $options

    $run_dir/runWorkflow.py -m local -j ${task.cpus}
	
    cp $run_dir/results/variants/genome.S1.vcf.gz $vcf
    cp $run_dir/results/variants/genome.S1.vcf.gz.tbi $tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS
    """

}
