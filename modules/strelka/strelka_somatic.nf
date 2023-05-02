process STRELKA_SOMATIC {

	tag "${meta.patient_id}|${meta.sample_id}"

        label 'long_parallel'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/STRELKA/", mode: 'copy'

	container 'quay.io/biocontainers/strelka:2.9.10--h9ee0642_1'

	input:
	tuple val(meta),path(normal),path(normal_bai),path(tumor),path(tumor_bai),path(indels),path(indels_tbi)
	tuple path(bed),path(bed_tbi)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf),path(tbi), emit: vcf
	tuple val(meta),path(vcf_indels),path(vcf_indels_tbi), emit: indels
	path("versions.yml"), emit: versions

	script:
	run_dir = "strelka_out"
	vcf = meta.patient_id + "_" + meta.sample_id + "-strelka_snps.vcf.gz"
	tbi = vcf + ".tbi"
	vcf_indels = meta.patient_id + "_" + meta.sample_id + "-strelka_indels.vcf.gz"
	vcf_indels_tbi = indels + ".tbi"

    """
    configureStrelkaSomaticWorkflow.py \
        --normalBam $normal \
        --tumorBam $tumor \
        --referenceFasta ${fasta} \
        --runDir $run_dir \
        --callRegions $bed \
        --indelCandidates $indels \
        --exome
		
    $run_dir/runWorkflow.py -m local -j ${task.cpus}

    # somatic.indels.vcf.gz  somatic.indels.vcf.gz.tbi  somatic.snvs.vcf.gz  somatic.snvs.vcf.gz.tbi
    cp $run_dir/results/variants/somatic.snvs.vcf.gz $vcf
    cp $run_dir/results/variants/somatic.snvs.vcf.gz.tbi $tbi
    cp $run_dir/results/variants/somatic.indels.vcf.gz $vcf_indels
    cp $run_dir/results/variants/somatic.indels.vcf.gz $vcf_indels_tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaGermlineWorkflow.py --version )
    END_VERSIONS
    """

}
