process GATK_APPLYBQSR {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/GATK", mode: 'copy'
    
    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    input:
    tuple val(meta),path(bam),path(bai),path(recal)
    path(intervals)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(recal_bam),path(recal_bai), emit: bam
    path("versions.yml"), emit: versions

    script:
    recal_bam = bam.getBaseName() + "-recal.cram"
    recal_bai = bam.getBaseName() + "-recal.cram.bai"

    """
    gatk  --java-options "-Xmx${task.memory.giga}g"  ApplyBQSR -R $fasta -I $bam -O $recal_bam -L $intervals -bqsr $recal \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
        --add-output-sam-program-record \
        --create-output-bam-md5 \
        --use-original-qualities -OBI true -ip $params.interval_padding

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
