//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowMain {

    //
    // Check and validate parameters
    //
    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        
        log.info header(workflow)

        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

    }

    public static String header(workflow) {
        def headr = ''
        def info_line = "IKMB Research Exome pipeline | version ${workflow.manifest.version}"
        headr = """
    ===============================================================================
    ${info_line}
    ===============================================================================
    """
        return headr
    }

    public static String help(workflow, params, log) {
        def command = "nextflow run ${workflow.manifest.name} --samples Samples.csv --assembly GRCh38 --kit xGen_v2 -profile diagnostic"
        def help_string = ''
        // Help message
        help_string = """

           Usage: nextflow run ikmb/exome-seq --assembly GRCh38 --kit xGen_v2 --samples Samples.csv
           This example will perform an exome analysis against the ALT-free hg38 assembly, assuming that exome reads were generated with
           the IDT xGen v2 kit and using DeepVariant with GLNexus.

           Required parameters:
           --samples                      A sample list in CSV format (see website for formatting hints)
           --assembly                     Name of the reference assembly to use
           --tools                        Comma-separated list of tools to run. 
           --kit                          Name of the exome kit (available options: xGen, xGen_custom, xGen_v2, Nextera, Pan_cancer)
           --email                        Email address to send reports to (enclosed in '')
           Optional parameters:
           --bwa2                         Use BWA2 instead of BWA1
           --phase                        Perform phasing of VCF files (only for multi-sample VCF files)
           --joint_calling                Perform joint calling of all samples (default: true)
           --skip_multiqc                 Don't attached MultiQC report to the email.
           --panel                        Gene panel to check coverage of (valid options: cardio_dilatative, cardio_hypertrophic, cardio_non_compaction, eoIBD_25kb, imm_eoIBD_full, breast_cancer)
           --all_panels                   Run all gene panels defined for this assembly (none if no panel is defined!)
           --panel_intervals              Run a custom gene panel in interval list format (must have a matching sequence dictionary!)
           --run_name                     A descriptive name for this pipeline run
           --cram                         Whether to output the alignments in CRAM format (default: bam)
           --interval_padding             Include this number of nt upstream and downstream around the exome targets (default: 10)
           --vep                          Perform variant annotation with VEP (requires substantial local configuration work!)
           Output:
           --outdir                       Local directory to which all output is written (default: results)
        """
        return help_string
    }

}
