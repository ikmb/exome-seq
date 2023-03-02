process CNVKIT_REFERENCE_TO_TARGETS {

    container 'quay.io/biocontainers/cnvkit:0.9.9--pyhdfd78af_0'

    input:
    path(cnn)

    output:
    path(bed), emit: targets

    script:

    bed = cnn.getBaseName() + ".bed"

    """
        reference2targets.py $cnn -o $bed
    """

}