process STRELKA {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/STRELKA", mode: 'copy'

	label 'strelka'

	input:
	tuple val(meta),path(bam),path(bai)
	tuple path(bed),path(bed_tbi)

	output:
	tuple val(meta),path(vcf), emit: vcf

	script:

	run_dir = "run_dir"
	vcf = meta.patient_id + "_" + meta.sample_id + "-strelka.vcf.gz"
        tbi = vcf + ".tbi"

	"""
		configureStrelkaGermlineWorkflow.py \
		--bam $bam \
		--referenceFasta ${params.fasta} \
		--runDir $run_dir \
		--callRegions $bed \
		--exome
		
		$run_dir/runWorkflow.py -m local -j ${task.cpus}

		cp $run_dir/results/variants/genome.S1.vcf.gz $vcf
                cp $run_dir/results/variants/genome.S1.vcf.gz.tbi $tbi

	"""

}

process STRELKA_JOINTCALLING {

        publishDir "${params.outdir}/MergedCallset/STRELKA_JOINT_CALLING", mode: 'copy'

        label 'strelka'

        input:
        tuple path(bams),path(bais)
        tuple path(bed),path(bed_tbi)

        output:
        tuple path(vcf_merged),path(vcf_merged_tbi), emit: vcf

        script:
        run_dir = "run_dir"
        vcf_merged = "strelka-joint_calling.vcf.gz"
	vcf_merged_tbi = vcf_merged + ".tbi"

        """
                configureStrelkaGermlineWorkflow.py \
                --bam ${bams.join(' --bam ')} \
                --referenceFasta ${params.fasta} \
                --runDir $run_dir \
                --callRegions $bed \
                --exome

                $run_dir/runWorkflow.py -m local -j ${task.cpus}

                cp $run_dir/results/variants/variants.vcf.gz $vcf_merged
		cp $run_dir/results/variants/variants.vcf.gz.tbi $vcf_merged_tbi

        """

}


