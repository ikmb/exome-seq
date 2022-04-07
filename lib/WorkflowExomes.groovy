//
// This file holds several functions specific to the workflow/esga.nf in the nf-core/esga pipeline
//

class WorkflowExomes {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {


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

    }

}
