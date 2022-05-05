process cnvkit_ref_decompress {
	
	input:
	path(cnn_gz)

	output:
	path(cnn)

	script:
	cnn = cnn_gz.getBaseName()

	"""
		gunzip -c $cnn_gz > $cnn
	"""

}

process cnvkit_ref_to_targets {


	label 'cnvkit'

	//publishDir "${params.outdir}/CnvKit/Ref", mode: 'copy'

	input:
	path(cnn)

	output:
	tuple path(targets),path(antitargets)

	script:
	targets = "cnvkit.target.bed"
	antitargets = "cnvkit.antitarget.bed"

	"""
		reference2targets.py $cnn -o cnvkit
	"""

}

// get per sample coverage for targets and antitargets
process cnvkit_coverage {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(meta), path(bam),path(bai)
	tuple path(targets),path(antitargets)

	output:
	tuple path(cnn),path(cnn_anti)
	tuple val(meta),path(cnn),path(cnn_anti)

	script:
	cnn = bam.getBaseName() + ".targetcoverage.cnn"
	cnn_anti = bam.getBaseName() + ".antitargetcoverage.cnn"

	"""
		cnvkit.py coverage $bam $targets -o $cnn
		cnvkit.py coverage $bam $antitargets -o $cnn_anti
	"""

}

// correct biases and stuff
process cnvkit_process {

	label 'cnvkit'
	input:
	tuple val(meta),path(cnn),path(cnn_anti)
	path(ref)

	output:
	tuple val(meta),path(cnr),path(cns)

	script:
	cns = cnn.getBaseName() + ".cns"
	cnr = cnn.getBaseName() + ".cnr"
	
	"""
		cnvkit.py fix $cnn $cnn_anti $ref -o $cnr
		cnvkit.py segment -p ${task.cpus} -m ${params.cnvkit_mode} $cnr -o $cns
	"""

}

// Segmetrics
process cnvkit_segmetrics {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(meta),path(cnr),path(cns)

	output:
	tuple val(meta),path(cnr),path(seg_cns)

	script:
	seg_cns = cns.getBaseName() + ".segmetrics.cns"

	"""
		cnvkit.py segmetrics -s $cns $cnr --ci --pi
	"""
}

// Attach confidence interfals
process cnvkit_call {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Processing", mode: 'copy'

	input:
	tuple val(meta),path(cnr),path(cns)

	output:
	tuple val(meta),path(cnr),path(call_cns)

	script:
	call_cns = cns.getBaseName() + ".call.cns"

	"""
		cnvkit.py call $cns -m threshold --filter ci
	"""

}

// Metrics per gene
process cnvkit_genemetrics {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Metrics", mode: 'copy'

	input:
	tuple val(meta), path(cnr),path(cns)

	output:
	path(metrics)

	script:

	metrics = cnr.getBaseName() + ".genemetrics.txt"

	"""
		cnvkit.py genemetrics -s $cns $cnr -t 0.2 -m 4 > $metrics
	"""

}

// find putative breaks
process cnvkit_breaks {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Metrics", mode: 'copy'

	input:
	tuple val(meta),path(cnr),path(cns)

	output:
	path(breaks)

	script:

	breaks = cnr.getBaseName() + ".breaks.txt"

	"""
		cnvkit.py breaks $cns $cnr > $breaks
	"""

}

// make useful output formats
process cnvkit_export {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit", mode: 'copy'

	input:
	tuple val(meta),path(cnr),path(call_cns)

	output:
	tuple path(bed),path(vcf)

	script:
	bed = call_cns.getBaseName() + ".bed"
	vcf = call_cns.getBaseName() + ".vcf"

	"""
		cnvkit.py export bed $call_cns -o $bed
		cnvkit.py export vcf $call_cns -i $meta.sample_id -o $vcf
	"""

}

process cnvkit_plots {

	label 'cnvkit'

	publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/CnvKit/Plots", mode: 'copy'

	input:
	tuple val(meta),path(cnr),path(call_cns)

	output:
	path(scatter)
	path(diagram)

	script:
	scatter = call_cns.getBaseName() + ".scatter.pdf"
	diagram = call_cns.getBaseName() + ".diagram.pdf"

	"""
		cnvkit.py scatter --y-min -4 --y-max 4 -o $scatter -s $call_cns $cnr
		cnvkit.py diagram -o $diagram -s $call_cns $cnr
	"""
}

