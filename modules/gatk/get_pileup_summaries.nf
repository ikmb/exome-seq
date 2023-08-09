process GATK_GET_PILEUP_SUMMARIES {

    tag "${meta.patient_id}|${meta.sample_id}"

    publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/MUTECT2/raw", mode: 'copy'

    container 'quay.io/biocontainers/gatk4:4.3.0.0--py36hdfd78af_0'

    label 'medium_serial'

    input:
    tuple val(meta),path(bam),path(bai)
    path(intervals)
    tuple path(fasta),path(fai),path(dict)

    output:
    tuple val(meta),path(stable), emit: table
    path("versions.yml"), emit: versions

    script:
    stable = bam.getBaseName() + "-pileup_summaries.table"

    """
    gatk GetPileupSummaries \
        -I $bam \
        -V $params.gnomad_af_vcf \
        -L $intervals \
        -R $fasta \
        -O $stable

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    
}
