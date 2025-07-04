process GATK_GENOTYPEGVCFS {

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'long_serial'

    input:
    tuple path(gvcf),path(tbi)
    path(intervals)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple path(vcf),path(tbi), emit: vcf
    path("versions.yml"), emit: versions

    script:
    vcf = gvcf.getSimpleName() + "-genotyped.vcf.gz"
    tbi = vcf + ".tbi"

    """
    gatk  --java-options "-Xmx${task.memory.giga}g" GenotypeGVCFs \
        -R $fasta \
        -V $gvcf -L $intervals -O $vcf \
        -G StandardAnnotation -G AS_StandardAnnotation \
        --allow-old-rms-mapping-quality-annotation-data \
        -ip $params.interval_padding -OVI true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
