process VALIDATE_AMPLICONBED {

        tag "${b}"

        input:
        path(b)

        output:
        path(ss), emit: bed

        script:
        ss = "amplicon_primers_validated.bed"

        """
                validate_ampliconbed.pl --infile $b
                cp $b $ss
        """
}
