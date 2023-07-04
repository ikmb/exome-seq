process GATK_VARIANTRECALIBRATOR {

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'short_serial'

    input:
    tuple val(meta),path(vcf),path(tbi)
    val(modus)

    output:
    tuple path(recal),path(recal_idx), emit: recal
    path(tranches), emit: tranches
    path("versions.yml"), emit: versions

    script:
    recal = modus + "-" + params.run_name + ".recal"
    recal_idx = recal + ".idx" 
    tranches = modus + "-" + params.run_name + ".tranches"

    def options = ""
    if (modus == "INDEL") {
        options = options.concat("-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR ")
        mills = params.genomes[ params.assembly ].mills        
        options = options.concat("--resource:mills,known=false,training=true,truth=true,prior=12 ${mills} ")
        axiom = params.genomes[ params.assembly ].axiom
        options = options.concat("--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom} ")
        dbsnp = params.genomes[ params.assembly ].dbsnp
        options = options.concat("--resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp} ")
    } else {
        options = options.concat("-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR ")
        options = options.concat("-max-gaussians 6 ")
        hapmap = params.genomes[ params.assembly ].hapmap
        options = options.concat("--resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap} ")
        omni = params.genomes[ params.assembly ].omni
        options = options.concat("--resource:omni,known=false,training=true,truth=true,prior=12 ${omni} ")
        g1k = params.genomes[ params.assembly ].g1k
        options = options.concat("--resource:1000G,known=false,training=true,truth=false,prior=10 ${g1k} ")
        dbsnp = params.genomes[ params.assembly ].dbsnp
        options = options.concat("--resource:dbsnp,known=true,training=false,truth=false,prior=7 ${dbsnp} ")
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
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

}
