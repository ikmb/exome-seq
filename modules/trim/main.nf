process trim {

        input:
        tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(center), val(date), path(fastqR1), path(fastqR2)

        output:
        tuple val(indivID), val(sampleID), val(libraryID), val(rgID), val(platform_unit), val(platform), val(platform_model), val(date), val(center), path(left),path(right)
        path(json)

        script:

        left = file(fastqR1).getBaseName() + "_trimmed.fastq.gz"
        right = file(fastqR2).getBaseName() + "_trimmed.fastq.gz"
        json = file(fastqR1).getBaseName() + ".fastp.json"
        html = file(fastqR1).getBaseName() + ".fastp.html"

        """
                fastp -c --in1 $fastqR1 --in2 $fastqR2 --out1 $left --out2 $right --detect_adapter_for_pe -w ${task.cpus} -j $json -h $html --length_required 35
        """
}

