process Strelka {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/Strelka", mode: 'copy'

	label 'strelka'

	input:
	tuple val(meta),path(bam)
	path(fasta)
	path(bed)

	output:
	tuple val(meta),path(vcf), emit: vcf
	script:

	run_dir = "run_dir"
	vcf = meta.patient_id + "_" + meta.sample_id + ".strelka." + params.run_name + ".vcf.gz"
	meta.caller = "Strelka"
	"""
		configureStrelkaGermlineWorkflow.py \
		--bam $bam \
		--referenceFasta $fasta \
		--runDir $run_dir \
		--callRegions $bed \
		--exome
		
		$run_dir/runWorkflow.py -m local -j ${task.cpus}

		cp $run_dir/results/variants/genome*.vcf.gz $vcf

	"""

}

process Strelka_JointCalling {

        publishDir "${params.outdir}/Strelka", mode: 'copy'

        label 'strelka'

        input:
        path(bams)
        path(fasta)
        path(bed)

        output:
        path(vcf), emit: vcf

        script:
        run_dir = "run_dir"
        vcf = "strelka." + params.run_name + ".vcf.gz"

        """
                configureStrelkaGermlineWorkflow.py \
                --bam ${bams.join(' --bam ')} \
                --referenceFasta $fasta \
                --runDir $run_dir \
                --callRegions $bed \
                --exome

                $run_dir/runWorkflow.py -m local -j ${task.cpus}

                cp $run_dir/results/variants/variants*.vcf.gz $vcf

        """

}

