// a simple sex check looking at coverage of the SRY gene
process sex_check {

        input:
        path(bams)

        output:
        path(sex_check_yaml)

        script:
        sex_check_yaml = "sex_check_mqc.yaml"

        """
                parse_sry_coverage.pl --fasta ${params.fasta} --region ${params.sry_bed} > $sex_check_yaml
        """
}

