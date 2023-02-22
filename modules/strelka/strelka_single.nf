process STRELKA_SINGLE_SAMPLE {

	tag "${meta.patient_id}|${meta.sample_id}"

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/STRELKA", mode: 'copy'

	label 'strelka'

	input:
	tuple val(meta),path(bam),path(bai)
	tuple path(bed),path(bed_tbi)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf), emit: vcf

	script:

	run_dir = "run_dir"
	vcf = meta.patient_id + "_" + meta.sample_id + "-strelka.vcf.gz"
        tbi = vcf + ".tbi"

	"""
		configureStrelkaGermlineWorkflow.py \
		--bam $bam \
		--referenceFasta ${fasta} \
		--runDir $run_dir \
		--callRegions $bed \
		--exome
		
		$run_dir/runWorkflow.py -m local -j ${task.cpus}

		cp $run_dir/results/variants/genome.S1.vcf.gz $vcf
                cp $run_dir/results/variants/genome.S1.vcf.gz.tbi $tbi

	"""

}
