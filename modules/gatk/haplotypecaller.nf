process GATK_HAPLOTYPECALLER {

	tag "${meta.patient_id}|${meta.sample_id}"
	
	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK/HC", mode: 'copy'

        container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

        label 'medium_serial'

	input:
	tuple val(meta),path(b),path(bai)
	path(intervals)
	val(modus)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(vcf),path(tbi), emit: vcf
	tuple val(meta),path(bam_out),path(bai_out), optional: true, emit: bam
	path("versions.yml"), emit: versions

	script:
	def options = ""
	if (modus == "single") {
		vcf = b.getBaseName() + "-hc.vcf.gz"
		tbi = vcf + ".tbi"
		bam_out = b.getBaseName() + "-hc.bam"
		bai_out = b.getBaseName() + "-hc.bai"
		options = "--bam-output $bam_out -OBI true"
	} else {
		vcf = b.getBaseName() + "-hc.vcf.gz"
		tbi = vcf + ".tbi"
		options = "-ERC GVCF -G AS_StandardAnnotation"
	}
	dbsnp = params.genomes[params.assembly].dbsnp
	
	//  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90
    """
    gatk HaplotypeCaller --java-options "-Xmx${task.memory.giga}g" -R $fasta -I $b -L $intervals -O $vcf \
        $options \
        -G StandardAnnotation -G StandardHCAnnotation \
        -OVI true -ip ${params.interval_padding} -D ${dbsnp} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}

