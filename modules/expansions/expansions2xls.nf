process EXPANSIONS2XLSX {

        publishDir "${params.outdir}/${meta.patient_id}/${meta.sample_id}/ExpansionHunter", mode: 'copy'

        input:
        tuple val(meta),path(report)

        output:
        path(expansion_xls)

        script:

        expansion_xls = report.getBaseName() + ".xlsx"

        """
                expansions2xls.pl --infile $report --outfile $expansion_xls
        """

}

