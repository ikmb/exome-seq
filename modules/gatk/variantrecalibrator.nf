process GATK_VARIANTRECALIBRATOR {

	label 'gatk'

	input:
	tuple val(meta),path(vcf),path(tbi)
	val(modus)

	output:
	tuple path(recal),path(recal_idx), emit: recal
	path(tranches), emit: tranches

	script:
	recal = modus + "-" + params.run_name + ".recal"
	recal_idx = recal + ".idx" 
	tranches = modus + "-" + params.run_name + ".tranches"

	def options = ""
	if (modus == "INDEL") {
		options = options.concat("-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP ")
		options = options.concat("--resource:mills,known=false,training=true,truth=true,prior=12 ${params.mills} ")
		options = options.concat("--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${params.axiom} ")
		options = options.concat("--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${params.dbsnp} ")
		options = options.concat("-max-gaussians 4 ")		
	} else {
		options = options.concat("-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP ")
		options = options.concat("-max-gaussians 6 ")
		options = options.concat("--resource:hapmap,known=false,training=true,truth=true,prior=15 ${params.hapmap} ")
		options = options.concat("--resource:omni,known=false,training=true,truth=true,prior=12 ${params.omni} ")
		options = options.concat("--resource:1000G,known=false,training=true,truth=false,prior=10 ${params.g1k} ")
		options = options.concat("--resource:dbsnp,known=true,training=false,truth=false,prior=7 ${params.dbsnp} ")
	}

	"""
		gatk VariantRecalibrator \
			--trust-all-polymorphic \
			-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
			$options \
			-mode $modus \
			-O $recal \
			--tranches-file $tranches -V $vcf

		gatk IndexFeatureFile -I $recal
	"""

}
