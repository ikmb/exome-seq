process CNVNATOR_EXTRACT {

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CNVnator", mode: 'copy'

	container 'quay.io/biocontainers/cnvnator:0.4.1--py310h0cf855d_6'

	input:
	tuple val(meta),path(bam),path(bai)
	tuple path(fasta),path(fai),path(dict)

	output:
	tuple val(meta),path(root), emit: root
	tuple val(meta),path(cnv_calls), emit: calls

	script:
	base_name = bam.getBaseName()
	root = bam.getBaseName() + ".root"
	cnv_calls = bam.getBaseName() + "-cnvnator.calls"

	"""
	cnvnator -root $root -tree $bam
	cnvnator -root $root -stat 1000 -fasta $fasta
	cnvnator -root $root -partition 1000
	cnvnator -root $root -call 1000 > $cnv_calls
	cnvnator2VCF.pl -prefix $base_name -reference GRCh38 -fasta $fasta $cnv_calls
	"""
}
