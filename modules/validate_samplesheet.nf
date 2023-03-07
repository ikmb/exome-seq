process VALIDATE_SAMPLESHEET {

        tag "${samplesheet}"

        input:
        path(samplesheet)

        output:
        path(ss), emit: csv

        script:
        ss = "Samples.validated.csv"

        """
                validate_samplesheet.pl --infile $samplesheet
                cp $samplesheet $ss
        """
}
