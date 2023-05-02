//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowExomes {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        
        genomeExistsError(params, log)
        kitExistsError(params, log)

        if (params.panel && params.all_panels || params.panel && params.panel_intervals || params.all_panels && params.panel_intervals) {
            log.info "The options for panel stats are mutually exclusive! Will use the highest ranked choice (panel > panel_intervals > all panels)"
        }

        if (!params.run_name) {
            log.info  "Must provide a run_name (--run_name)"
            System.exit(1)
        }
    
        if (!params.tools) {
            log.info "No analysis tools specified, performing only alignments and QC!"
        }

        if (params.build_references) {
            log.info "Building references...!"
        }

    }

    private static void genomeExistsError(params, log) {
        if (params.assembly && !params.genomes.containsKey(params.assembly)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome assembly '${params.assembly}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available assemblies  are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

    private static void kitExistsError(params, log) {
        if (params.kit && !params.genomes[params.assembly].kits.containsKey(params.kit)) {
            log.error "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Kit '${params.kit}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available kits for the chosen assembly are:\n" +
                "  ${params.genomes[params.assembly].kits.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            System.exit(1)
        }
    }

}
