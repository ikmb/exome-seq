process MANTA_TUMOR {

	tag "${meta.patient_id}|${meta.sample_id}"

	label 'medium_parallel'

	container 'quay.io/biocontainers/manta:1.6.0--h9ee0642_1'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MANTA", mode: 'copy'

	input:
	tuple val(meta),path(bam),path(bai)
	tuple path(bed_gz),path(bed_gz_tbi)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(sv),path(sv_tbi), emit: tumor_sv
	tuple val(meta),path(sv_can),path(sv_can_tbi), emit: candidate_sv
	tuple val(meta),path(indel),path(indel_tbi), emit: small_indels
	path("manta")
	path("versions.yml"), emit: versions

	script:
	sv = "${meta.patient_id}_${meta.sample_id}-tumorSV.vcf.gz"
	sv_tbi = sv + ".tbi"
	indel = "${meta.patient_id}_${meta.sample_id}-candidateSmallIndels.vcf.gz"
	indel_tbi = indel + ".tbi"
	sv_can = "${meta.patient_id}_${meta.sample_id}-candidateSV.vcf.gz"
	sv_can_tbi = sv_can + ".tbi"

    """
    configManta.py --tumorBam $bam --referenceFasta ${fasta} --runDir manta --callRegions $bed_gz --exome

    manta/runWorkflow.py -j ${task.cpus}

    cp manta/results/variants/tumorSV.vcf.gz $sv
    cp manta/results/variants/tumorSV.vcf.gz.tbi $sv_tbi
    cp manta/results/variants/candidateSmallIndels.vcf.gz $indel
    cp manta/results/variants/candidateSmallIndels.vcf.gz.tbi $indel_tbi
    cp manta/results/variants/candidateSV.vcf.gz $sv_can
    cp manta/results/variants/candidateSV.vcf.gz.tbi $sv_can_tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$( configManta.py --version )
    END_VERSIONS

    """

}

